#pragma once

#include <OpenColorIO/OpenColorIO.h>
#include "XYZ.h"

namespace OCIO = OCIO_NAMESPACE;

namespace ImageHandler {
    OCIO::ConstConfigRcPtr config = OCIO::Config::CreateFromFile("C:\\Users\\Charlie\\Libraries\\OCIOConfigs\\AgX-main\\config.ocio");
    OCIO::ConstProcessorRcPtr processor = config->getProcessor("Linear Rec.709", "AGX Base sRGB");
    auto compute = processor->getDefaultCPUProcessor();
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
        float post_lums = pow(lums, 1.0 / gamma);
        return in * post_lums / lums;
    }

    XYZ post_process_pixel(XYZ lums) {  //lums is raw light data via RGB
        lums = apply_gamma(lums, 1);
        float pixel[3] = { lums[0],lums[1],lums[2] };
        compute->applyRGB(pixel);
        lums = XYZ(pixel[0], pixel[1], pixel[2]);
        lums = apply_gamma(lums, 2.2);

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
