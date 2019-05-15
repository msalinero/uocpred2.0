import pandas as pd
import re

spanishPlayers = pd.Series(
    ['Casillas', 'Ramos', 'Pique', 'Puyol', 'Capdevila', 'Xabi Alonso', 'Busquets', 'Xavi Hernandez', 'Pedrito',
     'Iniesta', 'Villa'], index=[1, 15, 3, 5, 11, 14, 16, 8, 18, 6, 7])
print("Spanish Football Players: \n%s" % spanishPlayers)

kk = "------------MAALYRPGLRLNWHGLSPLGWPSCRSIQTLRVLSGDLGQLPTGIRDFVEHSARLCQPEGIHICDGTEAENTATLTLLEQQGLIRKLPKYNNCWLARTDPKDVARVESKTVIVTPSQRDTVPLPPGGARGQLGNWMSPADFQRAVDERFPGCMQGRTMYVLPFSMGP--------VGSPLSRIGVQLTDSAYVVASMRIMTRLGTPVLQALGDGDFVKCLHSVGQPLTGQGEPVSQWPCNPEKTLIGHVPDQREIISFGSGYGGNSLLGKKCFALRIASRLARDEGWLAEHMLILGITSPAGKKRYVAAAFPSACGKTNLAMMRPALPGWKVECVGDDIAWMRFDSEGRLRAINPENGFFGVAPGTSATTNPNAMATIQSNTIFTNVAETSDGGVYWEGIDQPLPPGVTVTSWLGKPWKPG-----MCGGEGVAQPPGLSTLMVEKLSPQPPTIF---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"

n = len(re.findall('-+', kk.strip('-')))
print(len(n))


