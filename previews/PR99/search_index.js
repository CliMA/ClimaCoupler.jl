var documenterSearchIndex = {"docs":
[{"location":"couplerstate/#Coupler-State","page":"Coupler State","title":"Coupler State","text":"","category":"section"},{"location":"couplerstate/","page":"Coupler State","title":"Coupler State","text":"The coupler provides a space to store information being passed between coupled model components at their boundaries. During this exchange, the coupler manages ancillary operations such as regridding, unit conversions, filtering, etc.","category":"page"},{"location":"couplerstate/","page":"Coupler State","title":"Coupler State","text":"The ClimaCoupler defines a type CouplerState for a container variable that holds information about the field boundary values that are being used to couple components. Components can use a coupler_put! operation to  export a set of field values to a CouplerState variable. A coupler_get operation is used to retrieve a set of field values from a CplState variable. During this exchange, the coupler manages ancillary operations such as  regridding, unit conversions, or filtering.","category":"page"},{"location":"couplerstate/#Coupler-Object-API","page":"Coupler State","title":"Coupler Object API","text":"","category":"section"},{"location":"couplerstate/","page":"Coupler State","title":"Coupler State","text":"    ClimaCoupler.CouplerState\n    ClimaCoupler.coupler_add_field!\n    ClimaCoupler.coupler_put!\n    ClimaCoupler.coupler_get","category":"page"},{"location":"couplerstate/#ClimaCoupler.CouplerState","page":"Coupler State","title":"ClimaCoupler.CouplerState","text":"CouplerState()\n\nType for holding coupler \"state\". This is the namespace through which coupled components communicate. Its role is to provide a level of indirection so that components remain modular and so that any data communication, interpolation, reindexing/unit conversions and filtering  etc... can be embeded in the intermdediate coupling layer.\n\nA field is exported by one component and imported by one or more other components.\n\n\n\n\n\n","category":"type"},{"location":"couplerstate/#ClimaCoupler.coupler_add_field!","page":"Coupler State","title":"ClimaCoupler.coupler_add_field!","text":"coupler_add_field!(\n        coupler::CouplerState,\n        fieldname::Symbol,\n        fieldvalue,\n    )\n\nAdd a field to the coupler that is accessible with key fieldname. \n\nArguments\n\ncoupler: coupler object the field is added to.\nfieldname: key to access the field in the coupler.\nfieldvalue: data array of field values.\n\n\n\n\n\n","category":"function"},{"location":"couplerstate/#ClimaCoupler.coupler_put!","page":"Coupler State","title":"ClimaCoupler.coupler_put!","text":"coupler_put!(coupler::CouplerState, fieldname::Symbol, fieldvalue)\n\nSets coupler field fieldname to fieldvalue.\n\n\n\n\n\n","category":"function"},{"location":"couplerstate/#ClimaCoupler.coupler_get","page":"Coupler State","title":"ClimaCoupler.coupler_get","text":"coupler_get(coupler::CouplerState, fieldname::Symbol [, target_sim::AbstractSimulation])\n\nRetrieve data array corresponding to fieldname.\n\nIf a target_sim is passed, the field is remapped to that simulation's boundary space.\n\n\n\n\n\n","category":"function"},{"location":"#ClimaCoupler.jl","page":"Home","title":"ClimaCoupler.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Coupling CliMA Models","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = ClimaCoupler","category":"page"},{"location":"","page":"Home","title":"Home","text":"ClimaCoupler.jl provides a means to couple CliMA  model components. It is designed to provide a flexible way to map boundary fluxes of quantities, like moisture and heat, that leave one component model (for example the atmosphere) to boundary fluxes of another component model (for example the ocean model). Functionality includes:","category":"page"},{"location":"","page":"Home","title":"Home","text":"coupled system time stepping control that integrates fluxes in time for sharing between components with differing time steps and/or time stepping schemes.\nsupport for mapping import and export boundary information between components so that fluxes of properties transferred between components are conserved.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The ClimaCoupler supports coupling components that are all within the same process or coupling components (using MPI) that are running on different processes.","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Coupler Scheme)\nClimaCoupler.jl allows for independent development of interchangeable component models.","category":"page"},{"location":"","page":"Home","title":"Home","text":"    ClimaCoupler","category":"page"},{"location":"#ClimaCoupler.ClimaCoupler","page":"Home","title":"ClimaCoupler.ClimaCoupler","text":"ClimaCoupler\n\nCoupling module sufficient for initial atmos-ocean-land coupled simulation.\n\n\n\n\n\n","category":"module"},{"location":"timestepping/#Coupled-Simulations-and-Timestepping","page":"Coupled Simulations & Timestepping","title":"Coupled Simulations & Timestepping","text":"","category":"section"},{"location":"timestepping/","page":"Coupled Simulations & Timestepping","title":"Coupled Simulations & Timestepping","text":"ClimaCoupler.jl organizes coupled models and their execution via the CoupledSimulation interface. An implementation of a CoupledSimulation, like any CliMA simulation, must  implement a ClimaCoupler.step! method, which in this context specifies the coupling scheme details (e.g. explicit, leap-frog, concurrent, etc...).","category":"page"},{"location":"timestepping/","page":"Coupled Simulations & Timestepping","title":"Coupled Simulations & Timestepping","text":"Component models being coupled must provide three methods for use in a  CoupledSimulation's ClimaCoupler.step! method:","category":"page"},{"location":"timestepping/","page":"Coupled Simulations & Timestepping","title":"Coupled Simulations & Timestepping","text":"step!: advances the component model in time\ncoupler_push!: prepares and puts coupled fields from the model   into the coupler via coupler_put! calls\ncoupler_pull!: gets coupled fields from the coupler   via coupler_get! calls and prepares them to be ingested by the model.","category":"page"},{"location":"timestepping/","page":"Coupled Simulations & Timestepping","title":"Coupled Simulations & Timestepping","text":"ClimaCoupler.CoupledSimulation\nClimaCoupler.run!\nClimaCoupler.step!\nClimaCoupler.coupler_push!\nClimaCoupler.coupler_pull!","category":"page"},{"location":"timestepping/#ClimaCoupler.run!","page":"Coupled Simulations & Timestepping","title":"ClimaCoupler.run!","text":"run!(::CoupledSimulation)\n\nA simple outer timestepping loop for coupled system runs.\n\nThis will be formalized when the run! functionality for component models is implemented so to have a consistent interface.\n\n\n\n\n\n","category":"function"},{"location":"timestepping/#ClimaCoupler.step!","page":"Coupled Simulations & Timestepping","title":"ClimaCoupler.step!","text":"step!(sim, dt)\n\nAdvances a simulation by dt.\n\nNote that dt is not necessarily the simulation's timestep length; a simuation could take several shorter steps that total to dt.\n\n\n\n\n\n","category":"function"},{"location":"timestepping/#ClimaCoupler.coupler_push!","page":"Coupled Simulations & Timestepping","title":"ClimaCoupler.coupler_push!","text":"coupler_push!(coupler::CouplerState, model)\n\nUpdate coupler with fields retrieved from the coupler.\n\ncoupler_push! is an adapter function to be implemented for each model component using the coupler. It should send coupling fields via coupler_put! calls and perform any operations on these fields to prepare them for the coupler.\n\n\n\n\n\n","category":"function"},{"location":"timestepping/#ClimaCoupler.coupler_pull!","page":"Coupled Simulations & Timestepping","title":"ClimaCoupler.coupler_pull!","text":"coupler_pull!(model, coupler::CouplerState)\n\nUpdate model with fields retrieved from the coupler.\n\ncoupler_pull! is an adapter function to be implemented for each model component using the coupler. It should get coupling fields via coupler_get calls and perform any operations on these fields to prepare them for use in the component model.\n\n\n\n\n\n","category":"function"}]
}