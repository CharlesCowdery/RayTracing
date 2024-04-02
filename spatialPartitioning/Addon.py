

bl_info = {
    "name": "SpatialRay Export Prepper",
    "description":"Preps scenes for export to SpatialRay"
}

import bpy

class PrepOBJ(bpy.types.Operator):
    bl_idname = "spatialray.prep"
    bl_label = "preps for export to SpaitalRay"
    bl_options = {"REGISTER"}
    def execute(self,context):
        for material in bpy.data.materials:
            if material.use_nodes:
                for node in material.node_tree.nodes:
                    if node.type == 'BSDF_PRINCIPLED':
                        if node.inputs[12].default_value == 0.5:
                            node.inputs[12].default_value = 0.5000001
                            break
        return {"FINISHED"}
    
def create(self, context):
    layout = self.layout
    scene = context.scene
    render = scene.render

    layout.separator()
    layout.operator(PrepOBJ.bl_idname, text="Prep Export")



def register():
    bpy.utils.register_class(PrepOBJ)
    bpy.types.TOPBAR_MT_edit.append(create)
    
def unregister():
    bpy.utils.register_class(PrepOBJ)

# make your unregister and remove draw
if __name__ == "__main__":
    register()