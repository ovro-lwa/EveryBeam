// Beam model coefficients converted by convert_coeff.py.
// Conversion performed on 2011/09/14/08:23:16 UTC using:
//     convert_coeff.py element_beam_HBA.coeff DefaultCoeffHBA.cc default_hba

#include <complex>

// Center frequency, and frequency range for which the beam model coefficients
// are valid. The beam model is parameterized in terms of a normalized
// frequency f' in the range [-1.0, 1.0]. The appropriate conversion is:
//
//     f' = (f - center) / range
//
const double default_hba_freq_center = 180e6;
const double default_hba_freq_range = 60e6;

// Shape of the coefficient array: 2x5x5x2 (the size of the last dimension is
// implied, and always equal to 2).
//
const unsigned int default_hba_coeff_shape[3] = {2, 5, 5};

// The array of coefficients in row-major order ("C"-order).
//
const std::complex<double> default_hba_coeff[100] = {
    std::complex<double>(0.9989322499459223, 0.0003305895124867),
    std::complex<double>(1.0030546028872600, 0.0002157249025076),
    std::complex<double>(0.0003002209532403, 0.0007909077657054),
    std::complex<double>(0.0022051270911392, 0.0003834815341981),
    std::complex<double>(-0.0003856663268042, 0.0008435910525861),
    std::complex<double>(0.0004887765294093, 0.0002777796480946),
    std::complex<double>(-0.0000699366665322, 0.0005136144371953),
    std::complex<double>(0.0001520602842105, 0.0001303481681886),
    std::complex<double>(0.0000512381993616, 0.0001550785137302),
    std::complex<double>(0.0000819244737818, 0.0000466470412396),
    std::complex<double>(0.0249658150445263, -0.0122024663463393),
    std::complex<double>(-0.0917825091832822, -0.0062606338208358),
    std::complex<double>(-0.0083709499453879, -0.0289759752488368),
    std::complex<double>(-0.0689260153643395, -0.0111348626546314),
    std::complex<double>(0.0116296166994115, -0.0307342946951178),
    std::complex<double>(-0.0171249717275797, -0.0080642275561593),
    std::complex<double>(0.0012408055399100, -0.0191295543986957),
    std::complex<double>(-0.0051031652662961, -0.0037143632875100),
    std::complex<double>(-0.0022414352263751, -0.0060474723525871),
    std::complex<double>(-0.0024377933436567, -0.0012852163337395),
    std::complex<double>(-0.6730977722052307, 0.0940030437973656),
    std::complex<double>(0.3711597596859299, 0.0557089394867947),
    std::complex<double>(0.2119250520015808, 0.2155514942677135),
    std::complex<double>(0.6727380529527980, 0.0989550572104158),
    std::complex<double>(-0.0419944347289523, 0.2355624543349744),
    std::complex<double>(0.1917656461134636, 0.0732470381581913),
    std::complex<double>(0.0048918921441903, 0.1588912409502319),
    std::complex<double>(0.0575369727210951, 0.0344677222786687),
    std::complex<double>(0.0241014578366618, 0.0547046570516960),
    std::complex<double>(0.0219986510834463, 0.0112189146988984),
    std::complex<double>(0.0665319393516388, -0.1418009730472832),
    std::complex<double>(-0.7576728614553603, -0.0472040122949963),
    std::complex<double>(-0.1017024786435272, -0.3302620837788515),
    std::complex<double>(-0.5600906156274197, -0.0797555201430585),
    std::complex<double>(0.0889729243872774, -0.3406964719938829),
    std::complex<double>(-0.1342560801672904, -0.0515926960946038),
    std::complex<double>(-0.0149335262655201, -0.2084962323582034),
    std::complex<double>(-0.0327252678958813, -0.0172371907472848),
    std::complex<double>(-0.0362395089905272, -0.0661322227928722),
    std::complex<double>(-0.0141568558526096, -0.0042676979206835),
    std::complex<double>(0.1121669548152054, 0.0504713119323919),
    std::complex<double>(0.1882531376700409, 0.0088411256350159),
    std::complex<double>(0.0066968933526899, 0.1181452711088882),
    std::complex<double>(0.0981630367567397, 0.0129921405004959),
    std::complex<double>(-0.0347327225501659, 0.1186585563636635),
    std::complex<double>(0.0102831315790362, 0.0046275244914932),
    std::complex<double>(0.0070209144233666, 0.0689639468490938),
    std::complex<double>(-0.0020239346031291, -0.0025499069613344),
    std::complex<double>(0.0132702874173192, 0.0207916487187541),
    std::complex<double>(0.0004387107229914, -0.0017223838914815),
    std::complex<double>(-0.0004916757488397, 0.0000266213616248),
    std::complex<double>(0.0006516553273188, -0.0000433166563288),
    std::complex<double>(-0.0004357897643121, 0.0000320567996700),
    std::complex<double>(0.0005818285824826, -0.0001021069650381),
    std::complex<double>(-0.0001047488648808, -0.0000302146563592),
    std::complex<double>(0.0001593350153828, -0.0000879125663990),
    std::complex<double>(-0.0000141882506567, -0.0000941521783975),
    std::complex<double>(-0.0000004226298134, -0.0000245060763932),
    std::complex<double>(-0.0000177429496833, -0.0000561890408003),
    std::complex<double>(-0.0000018388829279, 0.0000032387726477),
    std::complex<double>(0.0162495046881796, -0.0010736997976255),
    std::complex<double>(-0.0175635905033026, 0.0012997068962173),
    std::complex<double>(0.0138897851110661, -0.0014876219938565),
    std::complex<double>(-0.0150211436594772, 0.0029712291209158),
    std::complex<double>(0.0031705620225488, 0.0004838463688512),
    std::complex<double>(-0.0034418973689263, 0.0024603729467258),
    std::complex<double>(0.0003028387544878, 0.0026905629457281),
    std::complex<double>(0.0006768121359769, 0.0005901486396051),
    std::complex<double>(0.0004634797107989, 0.0016976603895716),
    std::complex<double>(0.0003344773954073, -0.0001499932789294),
    std::complex<double>(-0.1492097398080444, 0.0123735410547393),
    std::complex<double>(0.1393121453502456, -0.0121117146246749),
    std::complex<double>(-0.1217628319418324, 0.0222643129255504),
    std::complex<double>(0.1108579917761457, -0.0262986164183475),
    std::complex<double>(-0.0273147374272124, 0.0098595182007132),
    std::complex<double>(0.0208992817013466, -0.0205929453727953),
    std::complex<double>(-0.0002152227668601, -0.0089220757225133),
    std::complex<double>(-0.0074792188817697, -0.0043562231368076),
    std::complex<double>(-0.0012019994038721, -0.0079939660050373),
    std::complex<double>(-0.0035807498769946, 0.0014801422733613),
    std::complex<double>(0.1567990061437258, -0.0143275575385193),
    std::complex<double>(-0.1043118778001582, 0.0106756004832779),
    std::complex<double>(0.1151024257152241, -0.0225518489392044),
    std::complex<double>(-0.0593437249231851, 0.0216080058910987),
    std::complex<double>(0.0142781186223020, -0.0057037138045721),
    std::complex<double>(0.0151043140114779, 0.0141435752121475),
    std::complex<double>(-0.0057143555179676, 0.0141142700941743),
    std::complex<double>(0.0251435557201315, -0.0005753615445942),
    std::complex<double>(0.0004475745352473, 0.0102135659618127),
    std::complex<double>(0.0090474375150397, -0.0032177128650026),
    std::complex<double>(-0.0459124372023251, 0.0044990718645418),
    std::complex<double>(0.0135433541303599, -0.0021789296923529),
    std::complex<double>(-0.0306136798186735, 0.0064963361606382),
    std::complex<double>(-0.0046440676338940, -0.0037281688158807),
    std::complex<double>(-0.0006372791846825, 0.0008894047150233),
    std::complex<double>(-0.0181611528840412, -0.0011106177431486),
    std::complex<double>(0.0032325387394458, -0.0048123509184894),
    std::complex<double>(-0.0136340313457176, 0.0021185000810664),
    std::complex<double>(0.0001287985092565, -0.0032079544559908),
    std::complex<double>(-0.0045503800737417, 0.0015366231416036)};
