from abaqus import *
from abaqusConstants import *


def create_material(material_type, name_model, Frequency_idx, DensMod=0):
    """
    Create and set the material properties based on the material type.

    :return: Material name, density, elastic properties, and plastic properties
    """
    materials = {
        0: {
            'name': 'Ti-6Al-4V',
            'density': 4.429e-09,
            'elastic': (104000.0, 0.35),
            'plastic': ((
                            915.302610315806, 0.0), (923.782305619989, 0.00252525252525253), (
                            931.858019059606, 0.00505050505050505), (939.567323242512,
                                                                     0.00757575757575758),
                        (946.942220019862, 0.0101010101010101), (
                            954.0102494094, 0.0126262626262626), (960.795330934008,
                                                                  0.0151515151515152),
                        (967.31841174419, 0.0176767676767677), (
                            973.597972776577, 0.0202020202020202), (979.650428962764,
                                                                    0.0227272727272727),
                        (985.49044924186, 0.0252525252525253), (
                            991.131215086419, 0.0277777777777778), (996.584631331203,
                                                                    0.0303030303030303),
                        (1001.86149960238, 0.0328282828282828), (
                            1006.97166213041, 0.0353535353535354), (1011.92412189499,
                                                                    0.0378787878787879),
                        (1016.72714369522, 0.0404040404040404), (
                            1021.38833972509, 0.0429292929292929), (1025.91474247018,
                                                                    0.0454545454545455),
                        (1030.31286715786, 0.047979797979798), (
                            1034.58876554518, 0.0505050505050505), (1038.74807248004,
                                                                    0.053030303030303),
                        (1042.796046399, 0.0555555555555556), (
                            1046.73760471052, 0.0580808080808081), (1050.57735484172,
                                                                    0.0606060606060606),
                        (1054.31962159085, 0.0631313131313131), (
                            1057.96847131807, 0.0656565656565657), (1061.52773341843,
                                                                    0.0681818181818182),
                        (1065.00101944912, 0.0707070707070707), (
                            1068.391740224, 0.0732323232323232), (1071.70312113992,
                                                                  0.0757575757575758),
                        (1074.93821595962, 0.0782828282828283), (
                            1078.09991924234, 0.0808080808080808), (1081.1909775861,
                                                                    0.0833333333333333),
                        (1084.21399982213, 0.0858585858585859), (
                            1087.1714662825, 0.0883838383838384), (1090.0657372458,
                                                                   0.0909090909090909),
                        (1092.89906065158, 0.0934343434343434), (
                            1095.67357916251, 0.095959595959596), (1098.39133664338,
                                                                   0.0984848484848485),
                        (1101.05428411698, 0.101010101010101), (
                            1103.66428525012, 0.103535353535354), (1106.22312141611,
                                                                   0.106060606060606),
                        (1108.73249637493, 0.108585858585859), (
                            1111.19404060737, 0.111111111111111), (1113.60931533547,
                                                                   0.113636363636364),
                        (1115.97981625771, 0.116161616161616), (
                            1118.30697702463, 0.118686868686869), (1120.59217247756,
                                                                   0.121212121212121),
                        (1122.83672167087, 0.123737373737374), (
                            1125.04189069589, 0.126262626262626), (1127.2088953231, 0.128787878787879),
                        (1129.33890347717, 0.131313131313131), (1131.43303755822,
                                                                0.133838383838384),
                        (1133.49237662138, 0.136363636363636), (
                            1135.51795842547, 0.138888888888889), (1137.51078136074,
                                                                   0.141414141414141),
                        (1139.47180626449, 0.143939393939394), (
                            1141.40195813295, 0.146464646464646), (1143.30212773663,
                                                                   0.148989898989899),
                        (1145.17317314614, 0.151515151515152), (
                            1147.01592117443, 0.154040404040404), (1148.83116874139,
                                                                   0.156565656565657),
                        (1150.61968416583, 0.159090909090909), (
                            1152.38220838968, 0.161616161616162), (1154.11945613888,
                                                                   0.164141414141414),
                        (1155.83211702487, 0.166666666666667), (
                            1157.52085659046, 0.169191919191919), (1159.18631730367,
                                                                   0.171717171717172),
                        (1160.82911950244, 0.174242424242424), (
                            1162.44986229345, 0.176767676767677), (1164.04912440754,
                                                                   0.179292929292929),
                        (1165.62746501445, 0.181818181818182), (
                            1167.18542449915, 0.184343434343434), (1168.72352520185,
                                                                   0.186868686868687),
                        (1170.24227212398, 0.189393939393939), (
                            1171.74215360172, 0.191919191919192), (1173.22364194918,
                                                                   0.194444444444444),
                        (1174.68719407255, 0.196969696969697), (
                            1176.13325205704, 0.19949494949495), (1177.56224372789, 0.202020202020202),
                        (1178.97458318687, 0.204545454545455), (1180.37067132552,
                                                                0.207070707070707),
                        (1181.75089631628, 0.20959595959596), (
                            1183.11563408273, 0.212121212121212), (1184.46524874985,
                                                                   0.214646464646465),
                        (1185.8000930754, 0.217171717171717), (
                            1187.12050886329, 0.21969696969697), (1188.42682735982, 0.222222222222222),
                        (1189.7193696336, 0.224747474747475), (1190.99844694, 0.227272727272727), (
                            1192.26436107073, 0.22979797979798), (1193.5174046894, 0.232323232323232),
                        (1194.75786165351, 0.234848484848485), (1195.98600732376,
                                                                0.237373737373737),
                        (1197.20210886099, 0.23989898989899), (
                            1198.40642551147, 0.242424242424242), (1199.599208881, 0.244949494949495),
                        (1200.78070319841, 0.247474747474747), (1201.95114556872, 0.25)),
        },
        1: {
            'name': 'VeroClear',
            'density': 1.18e-09,
            'elastic': (1800.0, 0.30),
            'plastic': ((46, 0.0), (60.26, 0.02123), (63, 0.0885)),
        },
        2: {
            'name': 'VeroClearV2',
            'density': 1.18e-09,
            'elastic': (1013.0, 0.30),
            'plastic': ((52.392080003635236, 0.0), (52.62214293261519, 0.0002478526579665283),
                        (52.85242085578619, 0.0004956435753898092),
                        (53.06717011384611, 0.0007433883244012465), (53.27103766074756, 0.0009910823481716277),
                        (53.48927024040639, 0.001238700754420341),
                        (53.70771041318303, 0.001486257549347128), (53.92026837100669, 0.0017337587750301947),
                        (54.12028812872805, 0.0019812110324372023),
                        (54.328752999029525, 0.0022285936379143473), (54.535109198693355, 0.0024759170401611463),
                        (54.7334461275554, 0.00272318710425741),
                        (54.93181990316168, 0.002970395907780622), (55.131549680059244, 0.0032175421787921862),
                        (55.33347532753538, 0.0034646251184649243),
                        (55.53580137308763, 0.0037116465294720485), (55.73657541348726, 0.0039586083693803065),
                        (55.93233550900012, 0.004205514085899281),
                        (56.111801684291045, 0.00445237484442073), (56.291498270161064, 0.0046991743628261515),
                        (56.47915479108518, 0.004945905040914932),
                        (56.66650632034724, 0.0051925750676564755), (56.84964900884673, 0.0054391883268867305),
                        (57.028320350427244, 0.005685745107805777),
                        (57.20672442598782, 0.005932241290325678), (57.38488850300798, 0.00617867687754714),
                        (57.56038191876008, 0.006425054298936618),
                        (57.73563982051105, 0.006671371180571373), (57.90085338219223, 0.006917637235351465),
                        (58.05951734625152, 0.007163849043313228),
                        (58.21811021218013, 0.0074100002390042954), (58.376612809880385, 0.007656090871234922),
                        (58.534267211247105, 0.007902121718069872),
                        (58.689029609342896, 0.00814809482692961), (58.843574844648685, 0.008394007587119508),
                        (58.998237950458886, 0.008639859697703638),
                        (59.14625785184338, 0.008885657862766731), (59.28873341992808, 0.009131401027292887),
                        (59.424339702406265, 0.009377090528963211),
                        (59.55703471157241, 0.009622722490268604), (59.65926770538218, 0.00986832413807609),
                        (59.78698436772207, 0.010113840274327789),
                        (59.91754218474119, 0.010359293280606918), (60.039172840305945, 0.010604694803826718),
                        (60.15824431548734, 0.0108500385873282),
                        (60.27581482864058, 0.011095323616052248), (60.388404011085235, 0.011340553355259887),
                        (60.49830794909966, 0.011585725567866516),
                        (60.60826351070665, 0.011830837581639858), (60.71840512807913, 0.012075889293365089),
                        (60.8268120464392, 0.012320882628617275),
                        (60.930872457040046, 0.012565820195136017), (61.02418996868621, 0.012810708336653494),
                        (61.09351535734324, 0.013055560161768953),
                        (61.17479783623213, 0.01330034021199545), (61.25659774603036, 0.013545059809546768),
                        (61.33894634968156, 0.013789718952902014),
                        (61.421398902925425, 0.014034318110416374), (61.50403718914521, 0.014278857230653613),
                        (61.58625866268128, 0.014523336937702885),
                        (61.66811387194533, 0.014767757210918414), (61.74899959408394, 0.015012118675015146),
                        (61.824879801561316, 0.015256425343436686),
                        (61.89278242532009, 0.015500680179301549), (61.94909987531167, 0.015744886773071218),
                        (62.00244420477957, 0.01598903665236015),
                        (62.05261271801786, 0.016233130046388378), (62.10288340490778, 0.01647716374833788),
                        (62.151189514167974, 0.016721139827520534),
                        (62.19763309281907, 0.016965058212260534), (62.23812921401258, 0.017208922964103396),
                        (62.276567129609354, 0.017452730272727812),
                        (62.31126447151665, 0.017696481827910504), (62.33062808587734, 0.017940189103005155),
                        (62.3642450414883, 0.018183822919594944),
                        (62.38888130667878, 0.018427406242480293), (62.41151521010559, 0.01867093221179983),
                        (62.42834826618808, 0.01891440460618546),
                        (62.43632332461008, 0.019157826472425588), (62.4431495975144, 0.019401190229092438),
                        (62.45125261659644, 0.019644493510610442),
                        (62.461484684174785, 0.019887735504423465), (62.4717701087711, 0.02013091828837274),
                        (62.478770000156075, 0.020374045187265738),
                        (62.48291192835333, 0.020617115807765193), (62.48590182525875, 0.020860128494555596),
                        (62.49898071179499, 0.02110308205124603),
                        (62.50898071179499, 0.021345980227290756), (62.51898071179499, 0.021588828035138508),
                        (62.528980711794986, 0.02183161904302297),
                        (62.538980711794984, 0.022074351519396726), (62.54898071179498, 0.022317035012775446),
                        (62.55898071179498, 0.022559665366519252),
                        (62.56898071179498, 0.022802231027779674), (62.578980711794976, 0.023044737516722327),
                        (62.588980711794974, 0.023287185240109018),
                        (62.59898071179497, 0.0235295938670938), (62.60898071179497, 0.023771935071796244),
                        (62.61898071179497, 0.02401423567159751),
                        (62.628980711794966, 0.02425647100904084), (62.638980711794964, 0.02449864570634068),
                        (62.64898071179496, 0.024740761332196506),
                        (62.65898071179496, 0.024982818812136195), (62.66898071179496, 0.02522482635807767),
                        (62.67898071179496, 0.025466775777864215),
                        (62.688980711794954, 0.025708664330542835), (62.69898071179495, 0.025950492001621944),
                        (62.70898071179495, 0.026192267291122465),
                        (62.71898071179495, 0.026433989850317918), (62.72898071179495, 0.02667565888249522),
                        (62.738980711794945, 0.026917268530667823),
                        (62.74898071179494, 0.027158829934079254), (62.75898071179494, 0.0274003242097713),
                        (62.76898071179494, 0.027641760283392075),
                        (62.77898071179494, 0.027883139507095325), (62.788980711794935, 0.02812446083526028),
                        (62.79898071179493, 0.028365724057681688),
                        (62.80898071179493, 0.02860693074228915), (62.81898071179493, 0.028848081068484216),
                        (62.82898071179493, 0.029089173393374787),
                        (62.838980711794925, 0.029330213191120044), (62.84898071179492, 0.029571207696204152),
                        (62.85898071179492, 0.02981212117325034),
                        (62.86898071179492, 0.030052985578255503), (62.87898071179492, 0.03029379261111982),
                        (62.888980711794915, 0.03053454170785133),
                        (62.89898071179491, 0.030775232895329643), (62.90898071179491, 0.0310158675303807),
                        (62.91898071179491, 0.031256449390274055),
                        (62.92898071179491, 0.03149698958585391), (62.938980711794905, 0.031737458028307075),
                        (62.9489807117949, 0.031977868586169966),
                        (62.9589807117949, 0.03221822141391281), (62.9689807117949, 0.032458516642556834),
                        (62.9789807117949, 0.032698754224437174),
                        (62.988980711794895, 0.03293893414516662), (62.99898071179489, 0.03317905704101427),
                        (63.00898071179489, 0.03341912353747038),
                        (63.01898071179489, 0.033659132458601665), (63.02898071179489, 0.03389908464407894),
                        (63.038980711794885, 0.034138980613532474),
                        (63.04898071179488, 0.034378810030412654), (63.05898071179488, 0.0346185836636754),
                        (63.06898071179488, 0.034858301335466954),
                        (63.07898071179488, 0.0350979643034776), (63.088980711794875, 0.03533757153699296),
                        (63.09898071179487, 0.03557712188584099),
                        (63.10898071179487, 0.03581661517338988), (63.11898071179487, 0.03605606139488379),
                        (63.12898071179487, 0.03629545334248768),
                        (63.138980711794865, 0.0365347759605336), (63.14898071179486, 0.036774041416826414),
                        (63.15898071179486, 0.03701324980348615),
                        (63.16898071179486, 0.03725240110537435), (63.17898071179486, 0.03749149530500634),
                        (63.188980711794855, 0.03773053241700708),
                        (63.19898071179485, 0.03796951250418373), (63.20898071179485, 0.038208435309351044),
                        (63.21898071179485, 0.03844730098141225),
                        (63.22898071179485, 0.03868610990494982), (63.238980711794845, 0.03892486177292312),
                        (63.24898071179484, 0.039163556879341574),
                        (63.25898071179484, 0.039402189127050416), (63.26898071179484, 0.03964076068230961),
                        (63.27898071179484, 0.03987928237258424),
                        (63.288980711794835, 0.04011774903363194), (63.29898071179483, 0.04035615114226636),
                        (63.30898071179483, 0.04059449921941159),
                        (63.31898071179483, 0.04083279321955706), (63.32898071179483, 0.04107102996883734),
                        (63.338980711794825, 0.04130920904939621),
                        (63.34898071179482, 0.041547331505997945), (63.35898071179482, 0.041785396817347956),
                        (63.36898071179482, 0.042023404479689636),
                        (63.37898071179482, 0.04226135534245666), (63.388980711794815, 0.04249924970362459),
                        (63.39898071179481, 0.042737087508255944),
                        (63.40898071179481, 0.04297487063645837), (63.41898071179481, 0.043212597323545474),
                        (63.42898071179481, 0.043450266550781685),
                        (63.438980711794805, 0.04368788022771848), (63.4489807117948, 0.04392543694981435),
                        (63.4589807117948, 0.04416293318516816),
                        (63.4689807117948, 0.04440036847196717), (63.4789807117948, 0.04463774738815393),
                        (63.488980711794795, 0.04487507087138136),
                        (63.49898071179479, 0.045112340086247094), (63.50898071179479, 0.045349553107210366),
                        (63.51898071179479, 0.04558670987038564),
                        (63.52898071179479, 0.045823809658800466), (63.538980711794785, 0.04606085281681445),
                        (63.54898071179478, 0.04629783987001003),
                        (63.55898071179478, 0.04653477084726018), (63.56898071179478, 0.04677164584704188),
                        (63.57898071179478, 0.047008464912803125),
                        (63.588980711794775, 0.04724522805761197), (63.59898071179477, 0.04748194282213496),
                        (63.60898071179477, 0.047718600532744695),
                        (63.61898071179477, 0.04795519611971971), (63.62898071179477, 0.04819173546604966),
                        (63.638980711794765, 0.04842821894928247),
                        (63.64898071179476, 0.048664641177429135), (63.65898071179476, 0.04890099659140054),
                        (63.66898071179476, 0.0491373038558235),
                        (63.67898071179476, 0.04937356035781021), (63.688980711794756, 0.04960976109504253),
                        (63.69898071179475, 0.04984590612933627),
                        (63.70898071179475, 0.05008199506298533), (63.71898071179475, 0.050318026809964445),
                        (63.72898071179475, 0.05055400292747529),
                        (63.738980711794746, 0.05078992346851227), (63.748980711794744, 0.05102577498626287),
                        (63.75898071179474, 0.051261576244377396),
                        (63.76898071179474, 0.05149732644216397), (63.77898071179474, 0.051733021531680926),
                        (63.788980711794736, 0.05196866107560198),
                        (63.798980711794734, 0.05220424391046309), (63.80898071179473, 0.05243976898851131),
                        (63.81898071179473, 0.05267525016830146),
                        (63.82898071179473, 0.05291066943861027), (63.838980711794726, 0.05314602882084306),
                        (63.848980711794724, 0.0533813331557758),
                        (63.85898071179472, 0.05361658644477796), (63.86898071179472, 0.05385178393251197),
                        (63.87898071179472, 0.054086920395599536),
                        (63.888980711794716, 0.054322001177121285), (63.898980711794714, 0.054557026687674565),
                        (63.90898071179471, 0.054791996623441974),
                        (63.91898071179471, 0.05502691011619329), (63.92898071179471, 0.05526176656833602),
                        (63.938980711794706, 0.05549656175289948),
                        (63.948980711794704, 0.05573130492931652), (63.9589807117947, 0.055965997474034426),
                        (63.9689807117947, 0.0562006351286643),
                        (63.9789807117947, 0.056435217933467854), (63.988980711794696, 0.05666974600637266),
                        (63.998980711794694, 0.056904212723061565),
                        (64.00898071179469, 0.057138628886227866), (64.0189807117947, 0.057373000101947305),
                        (64.0289807117947, 0.05760730840295418),
                        (64.03898071179471, 0.05784156095013601), (64.04898071179471, 0.058075757466403546),
                        (64.05898071179472, 0.05830989694491623),
                        (64.06898071179472, 0.05854398170548032), (64.07898071179473, 0.05877801220617823),
                        (64.08898071179473, 0.0590119889587665),
                        (64.09898071179474, 0.05924590445174187), (64.10898071179474, 0.05947977357334559),
                        (64.11898071179475, 0.05971359220890694),
                        (64.12898071179475, 0.05994734886855588), (64.13898071179476, 0.060181051109306646),
                        (64.14898071179476, 0.06041469869026646),
                        (64.15898071179477, 0.06064829165686625), (64.16898071179477, 0.06088183014416257),
                        (64.17898071179478, 0.06111534181538347),
                        (64.18898071179478, 0.06134876992361364), (64.19898071179479, 0.061582143729823985),
                        (64.2089807117948, 0.06181546315175325),
                        (64.2189807117948, 0.062048728643166656), (64.2289807117948, 0.062281940408104586),
                        (64.23898071179481, 0.06251508913845623),
                        (64.24898071179481, 0.06274818068566065), (64.25898071179482, 0.062981220683002),
                        (64.26898071179482, 0.06321419915481176),
                        (64.27898071179483, 0.06344713364925181), (64.28898071179484, 0.063680016709828),
                        (64.29898071179484, 0.06391284563746542),
                        (64.30898071179485, 0.0641456191897041), (64.31898071179485, 0.06437833400549847),
                        (64.32898071179486, 0.06461099329854489),
                        (64.33898071179486, 0.06484360383138353), (64.34898071179487, 0.06507616462006002),
                        (64.35898071179487, 0.06530867017837244),
                        (64.36898071179488, 0.065541116299433), (64.37898071179488, 0.06577350835185018),
                        (64.38898071179489, 0.06600584653182588),
                        (64.39898071179489, 0.06623813068360988), (64.4089807117949, 0.06647036044111311),
                        (64.4189807117949, 0.06670253479294326),
                        (64.4289807117949, 0.06693465509033428), (64.43898071179491, 0.06716672159523479),
                        (64.44898071179492, 0.06739872360850772),
                        (64.45898071179492, 0.06763068344036742), (64.46898071179493, 0.06786259295852279),
                        (64.47898071179493, 0.06809443349074072),
                        (64.48898071179494, 0.0683262211928226), (64.49898071179494, 0.0685579559275043),
                        (64.50898071179495, 0.06878963870119666),
                        (64.51898071179495, 0.06902126804466813), (64.52898071179496, 0.0692528437879682),
                        (64.53898071179496, 0.06948436576646767),
                        (64.54898071179497, 0.06971585567636027), (64.55898071179497, 0.06994726474103373),
                        (64.56898071179498, 0.07017862076304984),
                        (64.57898071179498, 0.07040992639798792), (64.58898071179499, 0.070641177989349),
                        (64.598980711795, 0.07087237605370175),
                        (64.608980711795, 0.07110352074115871), (64.618980711795, 0.07133461214979916),
                        (64.62898071179501, 0.07156565024922007),
                        (64.63898071179501, 0.07179663488393032), (64.64898071179502, 0.07202756613815856),
                        (64.65898071179502, 0.07225844406218487),
                        (64.66898071179503, 0.07248926870086693), (64.67898071179503, 0.07272004006494868),
                        (64.68898071179504, 0.07295075869312813),
                        (64.69898071179504, 0.07318142607390815)),
        },
    }

    material_data = materials.get(material_type)
    if DensMod != 0:
        material_data['name'] = 'Ti-6Al-4V_Mod'
    if Frequency_idx == 0:
        material_data['density'] = 1
    elif Frequency_idx == 20:
        material_data['density'] = DensMod

    name_material = material_data['name']
    mdb.models[name_model].Material(name=name_material)
    mdb.models[name_model].materials[name_material].Density(table=((material_data['density'],),))
    mdb.models[name_model].materials[name_material].Elastic(
        table=((material_data['elastic'][0], material_data['elastic'][1]),))
    mdb.models[name_model].materials[name_material].Plastic(scaleStress=None, table=material_data['plastic'])

    return name_material
