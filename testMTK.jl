
### Events ###
continuous_events = [
[t ~ 3600 * shift1] => [GLC ~ GLC_1, ACT ~ ACT_1, BM ~ BM_1], # first shift
[t ~ 3600 * shift2] => [GLC ~ GLC_2, ACT ~ ACT_2, BM ~ BM_2] # second shift
]
