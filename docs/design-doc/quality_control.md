# **Quality control**

1. conservation 
    - total energy (or potential temperature density if prognostic)
    - mass 
    - (momentum)
2. performance 
    - comparison to standalone runs
    - scalability on multiple threads, nodes
3. compliance checking 
     - units: automated by Unitful 
4. error analysis of coupling strategies
    - isolated sequential splitting (default), dribbling, forcing, hybrid
    - see [Chris Vogl's talk](https://agu.confex.com/agu/fm21/meetingapp.cgi/Paper/839500) for leading order error semi-discrete analysis, [Ubbiali et al., 21](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020MS002377) for theoretical stability and accuracy analysis   

# Plan
- implement 1-3 as part of test cases *(when needed)*
- 4: toy model implementation of coupling strategies + SCM implementations *(long term)*  