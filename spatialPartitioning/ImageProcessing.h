#pragma once

#include <OpenColorIO/OpenColorIO.h>
#include "XYZ.h"

namespace OCIO = OCIO_NAMESPACE;

namespace ImageHandler {
    OCIO::ConstConfigRcPtr config = OCIO::Config::CreateFromFile("C:\\Users\\Charlie\\Libraries\\OCIOConfigs\\AgX-main\\config.ocio");
    OCIO::ConstProcessorRcPtr processor = config->getProcessor("Linear Rec.709", "AGX Base sRGB");
    auto compute = processor->getDefaultCPUProcessor();
    XYZ rgb2hsv(XYZ in) //https://stackoverflow.com/a/6930407/13946283
    {
        XYZ         out;
        double      min, max, delta;

        min = in.X < in.Y ? in.X : in.Y;
        min = min < in.Z ? min : in.Z;

        max = in.X > in.Y ? in.X : in.Y;
        max = max > in.Z ? max : in.Z;

        out.Z = max;                                // v
        delta = max - min;
        if (delta < 0.00001)
        {
            out.Y = 0;
            out.X = 0; // undefined, maybe nan?
            return out;
        }
        if (max > 0.0) { // NOTE: if Max is == 0, this divide would cause a crash
            out.Y = (delta / max);                  // s
        }
        else {
            // if max is 0, then r = g = b = 0              
            // s = 0, h is undefined
            out.Y = 0.0;
            out.X = NAN;                            // its now undefined
            return out;
        }
        if (in.X >= max)                           // > is bogus, just keeps compilor happy
            out.X = (in.Y - in.Z) / delta;        // between yellow & magenta
        else
            if (in.Y >= max)
                out.X = 2.0 + (in.Z - in.X) / delta;  // between cyan & yellow
            else
                out.X = 4.0 + (in.X - in.Y) / delta;  // between magenta & cyan

        out.X *= 60.0;                              // degrees

        if (out.X < 0.0)
            out.X += 360.0;

        return out;
    }


    XYZ hsv2rgb(XYZ in)
    {
        double      hh, p, q, t, ff;
        long        i;
        XYZ         out;

        if (in.Y <= 0.0) {       // < is bogus, just shuts up warnings
            out.X = in.Z;
            out.Y = in.Z;
            out.Z = in.Z;
            return out;
        }
        hh = in.X;
        if (hh >= 360.0) hh = 0.0;
        hh /= 60.0;
        i = (long)hh;
        ff = hh - i;
        p = in.Z * (1.0 - in.Y);
        q = in.Z * (1.0 - (in.Y * ff));
        t = in.Z * (1.0 - (in.Y * (1.0 - ff)));

        switch (i) {
        case 0:
            out.X = in.Z;
            out.Y = t;
            out.Z = p;
            break;
        case 1:
            out.X = q;
            out.Y = in.Z;
            out.Z = p;
            break;
        case 2:
            out.X = p;
            out.Y = in.Z;
            out.Z = t;
            break;

        case 3:
            out.X = p;
            out.Y = q;
            out.Z = in.Z;
            break;
        case 4:
            out.X = t;
            out.Y = p;
            out.Z = in.Z;
            break;
        case 5:
        default:
            out.X = in.Z;
            out.Y = p;
            out.Z = q;
            break;
        }
        return out;
    }
    static float luminance(XYZ v)
    {
        return XYZ::dot(v, XYZ(0.2126, 0.7152, 0.0722));
    }

    static XYZ change_luminance(XYZ c_in, float l_out)
    {
        float l_in = luminance(c_in);
        return c_in * (l_out / l_in);
    }

    //float realistic_response(float f, float iso) // https://graphics-programming.org/resources/tonemapping/index.html

    static XYZ apply_gamma(XYZ in, float gamma) {
        float lums = luminance(in);
        if (lums == 0) return XYZ(0, 0, 0);
        float post_lums = pow(lums, 1.0 / gamma);
        return in * post_lums / lums;
    }
    static XYZ apply_exposure(XYZ in, float exposure) {
        return in * pow(2, exposure);
    }
    XYZ post_process_pixel(XYZ lums) {  //lums is raw light data via RGB
        //lums = apply_exposure(lums, 2.6);
        lums = apply_gamma(lums, 1);

        float pixel[3] = { lums[0],lums[1],lums[2] };
        compute->applyRGB(pixel);
        lums = XYZ(pixel[0], pixel[1], pixel[2]);
        lums = apply_gamma(lums, 2.2);

        XYZ hsv = rgb2hsv(lums);
        //hsv.Y = std::min(hsv.Y * 1.4f, 1.0f);
        lums = hsv2rgb(hsv);

        return XYZ::clamp(lums, 0, 1) * 255;//scale to 24 bit rgb
    }
    static void postProcessRaw(vector<vector<XYZ*>>* data) {
        for (vector<XYZ*>& data_row : *data) {
            for (XYZ*& pixel_ptr : data_row) {
                *pixel_ptr = ImageHandler::post_process_pixel(*pixel_ptr);
            }
        }
    }
};
