# CouplerMachine
Repo for Coupler Specific Shared Development

General approach in coupler is to aim for something minimal to start that can grow over time.

First pass is for a simple slab ocean treated as a separate component beneath an atmosphere with
fields exchanged via the coupler abstratction. 

For this there are two big parts

## adding a bondary-fields abstraction
     The bondary-fields abstraction for each component would be used to expose and share
     state that is needed by another component. 
     The approach we are thinking about for this would be an DG object on the same mesh as the 
     component it is attached to, but possibly only 2-d (xy) or 2.5-d(xy + a thin layer in z).
     A starting point for ideas on this can be the code in 
     
        [https://github.com/CliMA/ClimateMachine.jl/blob/83cb930177d27ad411611c1aaaad68aec38ce4ec/src/Ocean/SplitExplicit01/Communication.jl#L191]
        
     which shows some mapping between differnet dimension collections of DG fields.
       
     
## creating a driver and stepping recipe for advancing coupled components
 
