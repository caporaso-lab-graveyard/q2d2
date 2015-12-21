import os
import warnings
import yaml
import frontmatter
import shutil

from .util import data_dirs, get_install_path

_plugin_dir = "plugins"
_manifest_n = "manifest.yml"
_using_name = "links.yml"
_plugins = None

def get_local_plugin_dir():
    plugin_dir = os.path.join(data_dirs.user_data_dir, _plugin_dir)
    if not os.path.exists(plugin_dir):
        os.makedirs(plugin_dir)
        with open(os.path.join(plugin_dir, _using_name), mode='w') as fh:
            fh.write(yaml.dump({}))

    return plugin_dir

class Workflow(object):
    def __init__(self, fp, id, plugin=None):
        self.id = id
        self.plugin = plugin
        self.fp = fp

        self._workflow = frontmatter.load(fp)
        self.title = self._workflow['name']
        self.inputs = set(self._workflow['inputs'])
        self.outputs = set(self._workflow['outputs'])

    def create_workflow(self, study_id='.'):
        workflow_fp = os.path.join(study_id, "%s.md" % self.id)
        if not os.path.exists(workflow_fp):
            with open(workflow_fp, mode='w') as dest_fh:
                dest_fh.write(self._workflow.content)
        return workflow_fp

    def delete_workflow(self, study_id='.'):
        os.remove(os.path.join(study_id, "%s.md" % self.id))

class Plugin(object):
    def __init__(self, fp, linked=False, name=None):
        self.broken = True
        self.name = name
        self.fp = os.path.abspath(fp)
        self.linked = linked
        manifest_fp = os.path.join(fp, _manifest_n)
        if not os.path.isfile(manifest_fp):
            warnings.warn("The plugin at %r is missing its %r file." %
                          (fp, _manifest_n), RuntimeWarning)
        else:
            with open(manifest_fp, mode='r', encoding='utf8') as fh:
                self.manifest = yaml.load(fh)
            self.name = self.manifest['name']
            self.broken = False

            self.workflows = {}
            for name in self.manifest['workflows']:
                id = name.split('.')[0]
                self.workflows[id] = Workflow(os.path.join(self.fp, name), id,
                                              plugin=self)

def load_links():
    with open(os.path.join(get_local_plugin_dir(), _using_name), mode='r',
              encoding='utf8') as fh:
        links = yaml.safe_load(fh)
        if links is None:
            links = {}
    return links

def write_links(links):
    with open(os.path.join(get_local_plugin_dir(), _using_name), mode='w',
              encoding='utf8') as fh:
        fh.write(yaml.dump(links))

def initialize_plugins():
    global _plugins
    _plugins = {}
    for plugin in yield_plugins():
        if plugin.name in _plugins:
            raise RuntimeError("Duplicate plugins found: %r" % self.name)
        elif plugin.broken:
            warnings.warn("Broken plugin %r at %r" % (plugin.name, plugin.fp))
        else:
            _plugins[plugin.name] = plugin

def get_plugins():
    if _plugins is None:
        initialize_plugins()

    return _plugins

def get_workflows():
    if _plugins is None:
        initialize_plugins()

    for plugin in _plugins.values():
        for workflow in plugin.workflows.values():
            yield workflow

def install_plugin(directory, linked=False):
    new_plugin = Plugin(directory, linked=linked)
    if new_plugin.broken:
        return False
    local_plugins = get_local_plugin_dir()
    if linked:
        links = load_links()
        links[new_plugin.name] = new_plugin.fp
        write_links(links)
    else:
        _, dir_name = os.path.split(new_plugin.fp)
        shutil.copytree(directory, os.path.join(local_plugins, dir_name))
    return True

def uninstall_plugin(name):
    for plugin in reversed(list(yield_plugins())):
        if plugin.name == name:
            if plugin.linked:
                links = load_links()
                del links[name]
                write_links(links)
            else:
                shutil.rmtree(plugin.fp)
            return True
    return False

def yield_plugins():
    packaged_plugins = os.path.join(get_install_path(), _plugin_dir)
    plugins = [{'fp': os.path.join(packaged_plugins, plugin)}
               for plugin in next(os.walk(packaged_plugins))[1]]

    local_plugins = get_local_plugin_dir()
    plugins += [{'fp': os.path.join(local_plugins, plugin)}
                for plugin in next(os.walk(local_plugins))[1]]

    plugins += [{'name': name, 'fp': fp, 'linked': True}
                for name, fp in load_links().items()]

    for plugin_init in plugins:
        yield Plugin(**plugin_init)
