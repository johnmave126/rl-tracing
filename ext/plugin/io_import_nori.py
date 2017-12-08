
bl_info = {
    "name": "Import Nori scenes format",
    "author": "Shuhao Tan",
    "version": (0, 1),
    "blender": (2, 7, 9),
    "location": "File > Import > Nori importer (.xml)",
    "description": "Import Nori scenes format (.xml)",
    "warning": "",
    "wiki_url": "",
    "tracker_url": "",
    "category": "Import-Export"}

import bpy
import bmesh
import os
import math
import shutil
import re
from xml.etree import ElementTree as et
from mathutils import Matrix, Vector


# Main class exporter


class NoriReader:

    def verbose(self, text):
        print(text)

    def __init__(self, context, filepath):
        self.context = context
        self.filepath = filepath
        self.workingDir = os.path.dirname(self.filepath)

    def __str2float(self, value):
        return re.split('[,\s]+', value)

    ######################
    # tools private methods
    # (xml format)
    ######################
    def read(self):
        tree = et.parse(self.filepath)

        # Clear the scene
        for ob in bpy.context.scene.objects:
            ob.select = True
        bpy.ops.object.delete()

        camera = tree.find("./camera")
        bpy.ops.object.camera_add()
        camera_obj = bpy.context.selected_objects[0]
        nearClip = 1e-4
        farClip = 1e4
        fov = 30.0
        width = 1280
        height = 720
        for parameter in camera:
            name = parameter.get('name')
            if name == 'nearClip':
                nearClip = float(parameter.get('value'))
            elif name == 'farClip':
                farClip = float(parameter.get('value'))
            elif name == 'fov':
                fov = float(parameter.get('value'))
            elif name == 'width':
                width = int(parameter.get('value'))
            elif name == 'height':
                height = int(parameter.get('value'))
            elif name == 'toWorld':
                lookat = parameter.find("./lookat")
                target = tuple(float(x)
                               for x in self.__str2float(lookat.get('target')))
                origin = tuple(float(x)
                               for x in self.__str2float(lookat.get('origin')))
                up = tuple(float(x)
                           for x in self.__str2float(lookat.get('up')))

                vtarget = Vector(target)
                vorigin = Vector(origin)
                vup = Vector(up)

                vdir = vtarget - vorigin
                vdir.normalize()
                vup.normalize()
                vleft = vup.cross(vdir)
                vleft.normalize()
                newup = vdir.cross(vleft)
                newup.normalize()

                lookat_transform = Matrix([
                    (vleft.x, newup.x, vdir.x, vorigin.x),
                    (vleft.y, newup.y, vdir.y, vorigin.y),
                    (vleft.z, newup.z, vdir.z, vorigin.z),
                    (0, 0, 0, 1)
                ])
                euler = lookat_transform.to_euler()
                euler.z = euler.z
                camera_obj.location = vorigin
                camera_obj.rotation_euler = euler

        camera_obj.data.clip_end = farClip
        camera_obj.data.clip_start = nearClip
        camera_obj.data.lens_unit = 'FOV'
        camera_obj.data.angle = fov / 180.0 * math.pi
        bpy.context.scene.render.resolution_x = width
        bpy.context.scene.render.resolution_y = height

        camera_obj.select = False
        meshes = tree.findall("./mesh")

        # import meshes
        for mesh in meshes:
            obj = mesh.find("./string[@name='filename']")
            obj_path = obj.get("value")

            bpy.ops.import_scene.obj(
                filepath=os.path.join(self.workingDir, obj_path))
            me = bpy.context.selected_objects[0].data

            bm = bmesh.new()
            bm.from_mesh(me)
            bmesh.ops.scale(bm, vec=(1, -1, -1))

            transforms = mesh.find("./transform")
            for transform in transforms:
                transform_type = transform.tag
                values = [float(x) for x in re.split(
                    '[,\s]+', transform.get('value'))]
                if transform_type == 'matrix':
                    bmatrix = list(zip(*([iter(values)] * 4)))
                    bmesh.ops.transform(bm, matrix=Matrix(bmatrix))
                elif transform_type == 'scale':
                    bmesh.ops.scale(bm, vec=values)
                elif transform_type == 'translate':
                    bmesh.ops.translate(bm, vec=values)

            bm.to_mesh(me)
            bm.free()

            mat = bpy.data.materials.new(
                name=bpy.context.selected_objects[0].name)

            bsdf = mesh.find("./bsdf")
            bsdf_type = bsdf.get("type")
            if bsdf_type == 'diffuse':
                albedo = bsdf.find("./color[@name='albedo']")
                color = tuple(float(x)
                              for x in self.__str2float(albedo.get('value')))
                bpy.context.selected_objects[0].data.materials.append(mat)
                mat.diffuse_color = color

            emitter = mesh.find("./emitter")
            if emitter is not None:
                radiance = emitter.find("./color[@name='radiance']")
                color = max(float(x) for x in re.split(
                    '[,\s]+', radiance.get('value')))
                mat.emit = color / 1000

            bpy.context.selected_objects[0].select = False


######################
# blender code
######################
from bpy.props import StringProperty, IntProperty, BoolProperty
from bpy_extras.io_utils import ExportHelper


class NoriImporter(bpy.types.Operator, ExportHelper):
    """Export a blender scene into Nori scene format"""

    # add to menu
    bl_idname = "import.nori"
    bl_label = "Import Nori scene"

    # filtering file names
    filename_ext = ".xml"
    filter_glob = StringProperty(default="*.xml", options={'HIDDEN'})

    def execute(self, context):
        nori = NoriReader(context, self.filepath)
        nori.read()
        return {'FINISHED'}

    def invoke(self, context, event):
        # self.frame_start = context.scene.frame_start
        # self.frame_end = context.scene.frame_end

        wm = context.window_manager
        wm.fileselect_add(self)
        return {'RUNNING_MODAL'}


def menu_import(self, context):
    import os
    default_path = os.path.splitext(bpy.data.filepath)[0] + ".xml"
    self.layout.operator(NoriImporter.bl_idname,
                         text="Import Nori scenes...").filepath = default_path


# Register Nori exporter inside blender
def register():
    bpy.utils.register_module(__name__)
    bpy.types.INFO_MT_file_import.append(menu_import)


def unregister():
    bpy.utils.unregister_module(__name__)
    bpy.types.INFO_MT_file_import.remove(menu_import)

if __name__ == "__main__":
    register()
