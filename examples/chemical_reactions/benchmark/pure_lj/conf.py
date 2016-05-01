Npart              = 32724
rho                = 0.8442
L                  = pow(Npart/rho, 1.0/3.0)
box                = (L, L, L)
r_cutoff           = 2.5
skin               = 1.0
temperature        = 0.5
dt                 = 0.0025
epsilon            = 1.0
sigma              = 1.0

warmup_cutoff      = pow(2.0, 1.0/6.0)
warmup_nloops      = 100
warmup_isteps      = 200
total_warmup_steps = warmup_nloops * warmup_isteps
epsilon_start      = 0.1
epsilon_end        = 1.0
epsilon_delta      = (epsilon_end - epsilon_start) / warmup_nloops
capradius          = 0.6
equil_nloops       = 10
equil_isteps       = 10
eq_conf            = 'eq_conf.pck'