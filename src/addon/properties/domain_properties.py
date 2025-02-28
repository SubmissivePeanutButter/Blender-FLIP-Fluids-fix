# Blender FLIP Fluids Add-on
# Copyright (C) 2023 Ryan L. Guy
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

if "bpy" in locals():
    import importlib
    reloadable_modules = [
        'domain_render_properties',
        'domain_bake_properties',
        'domain_simulation_properties',
        'domain_cache_properties',
        'domain_particles_properties',
        'domain_surface_properties',
        'domain_whitewater_properties',
        'domain_world_properties',
        'domain_presets_properties',
        'domain_materials_properties',
        'domain_advanced_properties',
        'domain_debug_properties',
        'domain_stats_properties'
    ]
    for module_name in reloadable_modules:
        if module_name in locals():
            importlib.reload(locals()[module_name])

import bpy, os, math
from mathutils import Vector, Color
from bpy.props import (
        BoolProperty,
        EnumProperty,
        FloatProperty,
        IntProperty,
        PointerProperty,
        StringProperty
        )

from . import (
        domain_render_properties,
        domain_bake_properties,
        domain_simulation_properties,
        domain_cache_properties,
        domain_particles_properties,
        domain_surface_properties,
        domain_whitewater_properties,
        domain_world_properties,
        domain_presets_properties,
        domain_materials_properties,
        domain_advanced_properties,
        domain_stats_properties,
        domain_debug_properties,
        preset_properties
        )
from .. import types
from ..objects import flip_fluid_cache
from ..utils import version_compatibility_utils as vcu
from ..utils import api_workaround_utils


class FlipFluidDomainProperties(bpy.types.PropertyGroup):
    conv = vcu.convert_attribute_to_28

    render = PointerProperty(
            name="Domain Render Properties",
            description="",
            type=domain_render_properties.DomainRenderProperties,
            ); exec(conv("render"))
    bake = PointerProperty(
            name="Domain Bake Properties",
            description="",
            type=domain_bake_properties.DomainBakeProperties,
            ); exec(conv("bake"))
    simulation = PointerProperty(
            name="Domain Simulation Properties",
            description="",
            type=domain_simulation_properties.DomainSimulationProperties,
            ); exec(conv("simulation"))
    cache = PointerProperty(
            name="Domain Cache Properties",
            description="",
            type=domain_cache_properties.DomainCacheProperties,
            ); exec(conv("cache"))
    particles = PointerProperty(
            name="Domain Surface Properties",
            description="",
            type=domain_particles_properties.DomainParticlesProperties,
            ); exec(conv("particles"))
    surface = PointerProperty(
            name="Domain Surface Properties",
            description="",
            type=domain_surface_properties.DomainSurfaceProperties,
            ); exec(conv("surface"))
    whitewater = PointerProperty(
            name="Domain Whitewater Properties",
            description="",
            type=domain_whitewater_properties.DomainWhitewaterProperties,
            ); exec(conv("whitewater"))
    world = PointerProperty(
            name="Domain World Properties",
            description="",
            type=domain_world_properties.DomainWorldProperties,
            ); exec(conv("world"))
    presets = PointerProperty(
            name="Domain Presets Properties",
            description="",
            type=domain_presets_properties.DomainPresetsProperties,
            ); exec(conv("presets"))
    materials = PointerProperty(
            name="Domain Materials Properties",
            description="",
            type=domain_materials_properties.DomainMaterialsProperties,
            ); exec(conv("materials"))
    advanced = PointerProperty(
            name="Domain Advanced Properties",
            description="",
            type=domain_advanced_properties.DomainAdvancedProperties,
            ); exec(conv("advanced"))
    debug = PointerProperty(
            name="Domain Debug Properties",
            description="",
            type=domain_debug_properties.DomainDebugProperties,
            ); exec(conv("debug"))
    stats = PointerProperty(
            name="Domain Stats Properties",
            description="",
            type=domain_stats_properties.DomainStatsProperties,
            ); exec(conv("stats"))
    mesh_cache = PointerProperty(
            name="Domain Mesh Cache",
            description="",
            type=flip_fluid_cache.FlipFluidCache,
            ); exec(conv("mesh_cache"))
    property_registry = PointerProperty(
            name="Domain Property Registry",
            description="",
            type=preset_properties.PresetRegistry,
            ); exec(conv("property_registry"))

    domain_settings_tabbed_panel_view = EnumProperty(
            name="Domain Panel View",
            description="Select settings panel to display",
            items=types.domain_settings_panel,
            default='DOMAIN_SETTINGS_PANEL_SIMULATION',
            options={'HIDDEN'},
            ); exec(conv("domain_settings_tabbed_panel_view"))


    def initialize(self):
        self.simulation.initialize()
        self.cache.initialize()
        self.advanced.initialize()
        self.materials.initialize()
        self._initialize_cache()
        self._initialize_property_registry()
        self.presets.initialize()


    def dummy_initialize(self):
        self.simulation.initialize()
        self.materials.initialize()
        self._initialize_property_registry()


    def destroy(self):
        self._delete_cache()


    def get_property_from_path(self, path):
        elements = path.split(".")
        if elements[0] == "domain":
            elements.pop(0)

        prop = self
        for e in elements:
            if not hasattr(prop, e):
                return None
            prop = getattr(prop, e)
            if (isinstance(prop, Vector) or isinstance(prop, Color) or
                (hasattr(prop, "__iter__") and not isinstance(prop, str))):
                new_prop = []
                for x in prop:
                    new_prop.append(x)
                prop = new_prop
            elif hasattr(prop, "is_min_max_property"):
                prop = [prop.value_min, prop.value_max]
        return prop


    def set_property_from_path(self, path, value):
        elements = path.split(".")
        if elements[0] == "domain":
            elements.pop(0)
        identifier = elements.pop()

        prop_group = self
        for e in elements:
            if not hasattr(prop_group, e):
                return False
            prop_group = getattr(prop_group, e)

        if not hasattr(prop_group, identifier):
            return False

        try:
            prop = getattr(prop_group, identifier)
            if hasattr(prop, "is_min_max_property"):
                prop.value_min = value[0]
                prop.value_max = value[1]
            else:
                setattr(prop_group, identifier, value)
        except:
            return False
        return True


    def _initialize_cache(self):
        self.mesh_cache.initialize_cache_objects()


    def _delete_cache(self):
        domain_object = bpy.context.scene.flip_fluid.get_domain_object()
        if domain_object is None:
            return
        self.mesh_cache.delete_cache_objects()


    def refresh_property_registry(self):
        # Call to re-initialize properties. The object in the scene is not
        # guaranteed to contain all properties of the current version of the addon,
        # such as a domain imported from the asset browser that was created in an earlier
        # version.
        self._initialize_property_registry()


    def _initialize_property_registry(self):
        self.property_registry.clear()
        self.render.register_preset_properties(     self.property_registry, "domain.render")
        self.bake.register_preset_properties(       self.property_registry, "domain.bake")
        self.simulation.register_preset_properties( self.property_registry, "domain.simulation")
        self.cache.register_preset_properties(      self.property_registry, "domain.cache")
        self.particles.register_preset_properties(  self.property_registry, "domain.particles")
        self.surface.register_preset_properties(    self.property_registry, "domain.surface")
        self.whitewater.register_preset_properties( self.property_registry, "domain.whitewater")
        self.world.register_preset_properties(      self.property_registry, "domain.world")
        self.presets.register_preset_properties(    self.property_registry, "domain.presets")
        self.materials.register_preset_properties(  self.property_registry, "domain.materials")
        self.advanced.register_preset_properties(   self.property_registry, "domain.advanced")
        self.debug.register_preset_properties(      self.property_registry, "domain.debug")
        self.stats.register_preset_properties(      self.property_registry, "domain.stats")
        self._validate_property_registry()


    def _validate_property_registry(self):
        for p in self.property_registry.properties:
            path = p.path
            base, group, identifier = path.split('.', 2)
            if not hasattr(self, group):
                print("Property Registry Error: Unknown Property Group <" + 
                      group + ", " + path + ">")
                continue
            prop_group = getattr(self, group)
            if not hasattr(prop_group, identifier):
                print("Property Registry Error: Unknown Identifier <" + 
                      identifier + ", " + path + ">")
                continue


    def scene_update_post(self, scene):
        self.render.scene_update_post(scene)
        self.simulation.scene_update_post(scene)
        self.surface.scene_update_post(scene)
        self.world.scene_update_post(scene)
        self.debug.scene_update_post(scene)
        self.stats.scene_update_post(scene)
        self.materials.scene_update_post(scene)


    def frame_change_post(self, scene, depsgraph=None):
        api_workaround_utils.frame_change_post_apply_T71908_workaround(bpy.context, depsgraph)
        self.world.frame_change_post(scene)
        self.stats.frame_change_post(scene, depsgraph)


    def load_pre(self):
        pass


    def load_post(self):
        self.simulation.load_post()
        self.bake.load_post()
        self.cache.load_post()
        self.surface.load_post()
        self.stats.load_post()
        self.debug.load_post()
        self.presets.load_post()
        self.advanced.load_post()
        self.materials.load_post()
        self.mesh_cache.load_post()
        self._initialize_property_registry()

        api_workaround_utils.load_post_update_cycles_visibility_forward_compatibility_from_blender_3()


    def save_pre(self):
        self.debug.save_pre()


    def save_post(self):
        self.cache.save_post()


def scene_update_post(scene):
    dprops = scene.flip_fluid.get_domain_properties()
    if dprops is None:
        return
    if not scene.flip_fluid.is_domain_in_active_scene():
        return
    dprops.scene_update_post(scene)


def frame_change_post(scene, depsgraph=None):
    dprops = scene.flip_fluid.get_domain_properties()
    if dprops is None:
        return
    dprops.frame_change_post(scene, depsgraph)


def load_pre():
    dprops = bpy.context.scene.flip_fluid.get_domain_properties()
    if dprops is None:
        return
    dprops.load_pre()


def load_post():
    dprops = bpy.context.scene.flip_fluid.get_domain_properties()
    if dprops is None:
        return
    dprops.load_post()


def save_pre():
    dprops = bpy.context.scene.flip_fluid.get_domain_properties()
    if dprops is None:
        return
    dprops.save_pre()


def save_post():
    dprops = bpy.context.scene.flip_fluid.get_domain_properties()
    if dprops is None:
        return
    dprops.save_post()


def register():
    domain_render_properties.register()
    domain_bake_properties.register()
    domain_simulation_properties.register()
    domain_cache_properties.register()
    domain_particles_properties.register()
    domain_surface_properties.register()
    domain_whitewater_properties.register()
    domain_world_properties.register()
    domain_presets_properties.register()
    domain_materials_properties.register()
    domain_advanced_properties.register()
    domain_debug_properties.register()
    domain_stats_properties.register()
    bpy.utils.register_class(FlipFluidDomainProperties)


def unregister():
    domain_render_properties.unregister()
    domain_bake_properties.unregister()
    domain_simulation_properties.unregister()
    domain_cache_properties.unregister()
    domain_particles_properties.unregister()
    domain_surface_properties.unregister()
    domain_whitewater_properties.unregister()
    domain_world_properties.unregister()
    domain_presets_properties.unregister()
    domain_materials_properties.unregister()
    domain_advanced_properties.unregister()
    domain_debug_properties.unregister()
    domain_stats_properties.unregister()
    bpy.utils.unregister_class(FlipFluidDomainProperties)