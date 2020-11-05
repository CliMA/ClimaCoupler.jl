# Simple coupled column equations to play with simple coupled DG datatypes/code design
#
#   equations have one trivial source term and one trivial coupling term 
#   
#   equations couple 
#   
#   * atmos col temp over ocean ($T^{ao}$) and ocean col temp under atmos ($T^{oa}$)
#   
#   * atmos col temp over land ($T^{al}$) and land temp under atmos ($T^{la}$)
#   
#   $ \frac{d T^{ao} }{dt}= -\frac{1}{\tau^{aa}}\left( T^{ao}-T^{aref}(t) \right) -\frac{1}{\tau^{ao}}\left( T^{ao}-T^{oa}\right)$
#
#   $ \frac{d T^{oa} }{dt}= -\frac{1}{\tau^{oo}}\left( T^{oa}-T^{odeep} \right) +\frac{1}{\tau^{ao}}\left( T^{ao}-T^{oa}\right)$
#
#   $ \frac{d T^{al} }{dt}= -\frac{1}{\tau^{aa}}\left( T^{al}-T^{aref}(t) \right) -\frac{1}{\tau^{al}}\left( T^{al}-T^{la}\right)$
#
#   $ \frac{d T^{la} }{dt}= -\frac{1}{\tau^{ll}}\left( T^{la}-\alpha^l T^{aref}(t) \right) +\frac{1}{\tau^{al}}\left( T^{al}-T^{la}\right)$
#
#
#

# Some time constants (hour, day, year)
Hˢ=3600
Dˢ=86400
Yˢ=Dˢ*365

# A timestep
Δt=Hˢ

# Relaxation timescales
τᵃᵃ=6*Hˢ
τˡˡ=6*Hˢ
τᵒᵒ=5*Yˢ
τᵃᵒ=60*Dˢ
τᵃˡ=1*Dˢ

# How many steps
Nˢ=Int(30*Yˢ/Δt)

# Time dependent equilibrated "atmos" temperature
Tᵐⁱᵈ=10
Tᵃᵐᵖ=10
Tᵃʳᵉᶠ(t)=( Tᵐⁱᵈ + Tᵃᵐᵖ * sin(t/Yˢ*2*π  ) )
# Fixed backgroud "ocean" state
Tᵒᵈᵉᵉᵖ=4
αˡ=1.3

# Initial condition and allocate memory
Tᵃᵒ=ones(Nˢ)*Tᵃʳᵉᶠ(0)
Tᵒᵃ=ones(Nˢ)*Tᵃʳᵉᶠ(0)
Tᵃˡ=ones(Nˢ)*Tᵃʳᵉᶠ(0)
Tˡᵃ=ones(Nˢ)*Tᵃʳᵉᶠ(0)

for i=1:Nˢ-1
 tim=i*Δt
 s_AO=-1/τᵃᵃ*( Tᵃᵒ[i]-Tᵃʳᵉᶠ(tim)    )
 s_OA=-1/τᵒᵒ*( Tᵒᵃ[i]-Tᵒᵈᵉᵉᵖ        )
 s_LA=-1/τˡˡ*( Tˡᵃ[i]-Tᵃʳᵉᶠ(tim)*αˡ )
 s_AL=-1/τᵃᵃ*( Tᵃˡ[i]-Tᵃʳᵉᶠ(tim)    )

 f_A2O=-1/τᵃᵒ*(Tᵃᵒ[i]-Tᵒᵃ[i  ])
 f_A2L=-1/τᵃˡ*(Tᵃˡ[i]-Tˡᵃ[i  ])

 
 Tᵃᵒ[i+1]=Tᵃᵒ[i]+s_AO*Δt+f_A2O*Δt
 Tᵒᵃ[i+1]=Tᵒᵃ[i]+s_OA*Δt-f_A2O*Δt
 Tˡᵃ[i+1]=Tˡᵃ[i]+s_LA*Δt-f_A2L*Δt
 Tᵃˡ[i+1]=Tᵃˡ[i]+s_AL*Δt+f_A2L*Δt  
end

using Plots

plot(Tᵒᵃ,label="Tᵒᵃ")
