Plugins and Packages
====================

While CHARMMing was originally designed for CHARMM, in its third revision it was redesigned from the ground up to support a wide range of programs, such as Q-chem. Even better, while we have a few packages forthcoming that will give generic tasks for CHARMM, Q-chem, and more, you may want to create additional functionality surrounding specific features such as determining if an uploaded PDB exists in an internal database, and you can now use CHARMMing's build in plugin architecture to accomplish this. 

# Internal plugins

If you have a need to build a plugin that "blocks", i.e. something that happens syncronously and can stop the process of a specific workflow based on set events, you can create an internal plugin that listens to various signals that CHARMMing emits and hooks into them. It's also useful for asyncronous programs that need to access more internal CHARMMing specifics that's only possible by using Python. These plugins are detected automatically by CHARMMing and will show up in the admin interface as a loaded plugin, and if you have settings, they will appear with defaults if present. Note that while the plugins are loaded automatically and CHARMMing will catch any syntax or other breaking errors when loading the plugins, it's not possible to check each piece of functionality until it's run, so be wary of any code that could break CHARMMing during its execution. 

# External plugins

If you want to build a plugin that listens for signals from CHARMMing but doesn't necessarily need to interact with CHARMMing, i.e. a logger, you can create an external plugin. It uses JSON as its data type and uses either ZeroMQ or HTTP REST as its protocol. 