#include "include/fp.h"

// TODO should we modify fp_sqrt to take an out
// variable?
void 
fp_sqrt(fp_t * x)
{   
    fp_t tmp = *x;
    (void)gf27500_sqrt(x, &tmp);
}

bool
fp_is_square(const fp_t * a)
{
    return -1 != gf27500_legendre(a);
}

// TODO should we modify fp_inv to take an out
// variable?
void
fp_inv(fp_t * x)
{
    fp_t tmp = *x;
    (void)gf27500_invert(x, &tmp);
}

/*
CAUTION: HAZARDOUS MATERIALS

The following tables are constants used for torsion basis generation. All values have already been converted into Montgomery
form, and so if the internal API for the field changes, these values may also need to change

NQR_TABLE:
    is a table of 20 values all of which are not squares in GF(p^2) with modulus x^2 + 1 and prime p = 27*2**500 - 1
Z_NQR_TABLE:
    is a table of 20 values for which the value is a square and (value - 1) is not a square in the field GF(p^2) 
    with modulus x^2 + 1 and prime p = 27*2**500 - 1 
*/

const uint64_t NQR_TABLE[20][2][NWORDS_FIELD] = {
    {{0x9B9C9D16BC8714CD, 0xA8CCD3D5AA9B25BE, 0x5CB0B229AFFB031A, 0x1034E0FAD3879BFD, 0xADF66A044A5CB6D4, 0xCEFE95E3E38E0F77, 0xCCEAE5A848DF543D, 0x0027C7DE38F22242},
     {0x134DD7E8E79ED23D, 0xC2D74B1B59F759A3, 0xF83EE94DFB78CD22, 0x8ABFBBF2C794F2D5, 0x4A922D5252DE3F9E, 0xEACB297359A638CA, 0x047EA519584BDEAA, 0x00559B7C1FAA7C13}},
    {{0x7C89C343BFF4B3D6, 0xC481DB28411841BD, 0x0F6775A92F28A17D, 0x8C29716D84A01180, 0x2611C90B35B8DA6A, 0xDB3B400B9E242963, 0x120B8BDEF6C4BEBC, 0x014630EFC7FC3A9E},
     {0x959A8AC96BCC61B6, 0x0E3F7572877B9258, 0x7EC388FAA822B366, 0x8F222CFAA21A5785, 0x6ED1B70424239857, 0x8C7F659E663F870C, 0xA7EF3A71519C1483, 0x018457C415603C01}},
    {{0x332F42117D2812AD, 0x1CADB01013267AC3, 0x7DBE43C3326B45CE, 0x78A92B9EB2D09FF6, 0x88EE8D4CD764B796, 0x5CFDF295E14782B8, 0x1A5012C22BF55824, 0x005422E55591F171},
     {0xA5189DE329614F8E, 0x625B9755100A0E6D, 0x8E9AEAD617D4D154, 0xE2C1F14B0ED1A672, 0x056DE1B6353B14EB, 0x10846245CA402433, 0x6A8461732747E77B, 0x0084D8408BD7BD06}},
    {{0x9392283B6FC6A05A, 0x5B4ADBC32EFE9E37, 0xFD349C9A44110C1E, 0x4E6B993D4ABFA62F, 0xF737E20422B4EB2E, 0x3B1C158212530FDB, 0xFA0F69816870DF69, 0x00DE55AFCB6C1B04},
     {0x59F28166F5E665FA, 0xDF207B71A7E54832, 0x2CA5E023FA157715, 0xBBA416934752A94D, 0x9F95C1165C11EC95, 0x4CF1674A6886F77E, 0xC2B7D963ED43FE24, 0x004863D33A0F5FDF}},
    {{0xC6276224A3DA2EDF, 0x2FAD9F92FF9CF355, 0xD33F1C0ACA987699, 0x6D7C720ABECF9D1B, 0x16949D8976C370CA, 0x75841B6AB67B813B, 0x4619E398D45A9D67, 0x01486DEE05F2976E},
     {0x6521190117AFE2BF, 0x1B5780A522F0E23B, 0x7A96D1294F16C4EE, 0xC6E5859ECD28A088, 0x59D3C220A50905E0, 0xB960AB599A9244C7, 0xE0BB8CDE73C81A9A, 0x003E800053323EEA}},
    {{0xC31644D80089FBDE, 0x6E9DAD3D003A3C50, 0x7C61CB348803C0F5, 0x7D8A64F8C6E7636B, 0x346C754236F28E62, 0x251FE498A03856FA, 0x4010D0C16FF5C6AB, 0x019CE06CCF4EFFF6},
     {0x4EF4888CE08B3772, 0xC86D1C3FF28DB6A5, 0xF54BB25571876F38, 0x0D0C08D50AAD1C4C, 0x4D8DFE8B895FFB95, 0xFF7967883B55BADB, 0xC2D95A41D2F1F767, 0x0122EAE9085659FA}},
    {{0x2A1D4CB67AA37762, 0x80E99CB36527586F, 0xBB8669BCEFB1D157, 0x492F9971BAB573B6, 0x9AF6F705F52BCDFD, 0xED4CBAD9188F7DB1, 0x2016295D80C97305, 0x004B59FF51B2BA2C},
     {0xA44C0EE3BB40CCAE, 0x28A499A882B30641, 0xB4F0D76DE783687A, 0xD62B0313579D0318, 0xD53E453F74E870B5, 0xCA3248FB1FFF43D6, 0x2ECDA0EEF31614AD, 0x005CA71581293C32}},
    {{0xEFCDD17ABD2A95D5, 0x1B17F2D63E26834D, 0x086E646FA3AA46BC, 0x66A8D4E420BF9482, 0x002D573201EDEFC8, 0x0878DCEA1078198D, 0x75A1C0927A965777, 0x00A7F6B982E70093},
     {0x553C78373B0FD804, 0xC377ED8B4350E9BA, 0x93E95E32953C4691, 0xA42EF53B45D48552, 0x9F7880D2D772CDC4, 0x325DE75125A2698C, 0x4FDF8C30EA357ACF, 0x0142DD6D316174A7}},
    {{0xD2E5B8C0612EE65B, 0x15EBD106028A9A15, 0x80FDD2EF844000FA, 0x7A4E65DEA0E5EA4E, 0xF69C46D456918012, 0x2CEC358377C4283E, 0x28AFBF38BA8447DD, 0x00435A5EAC5F510D},
     {0x4E81365E9618C9FD, 0x5DA9ED490A66D76D, 0x30B18E2880A6D625, 0x7524DF3A9EF3454D, 0xCBB598CAABBE09B2, 0x13DB46D376CEE614, 0x0DF12145B2534357, 0x00CEE0232D9F2AA8}},
    {{0xD00848FFD46793CD, 0x97F25A86F93A6F78, 0x6E3D089ACB4988FE, 0x7232AF2D89C5A305, 0xE1FF2DEC79B56EE9, 0xC9F135492C70D78E, 0x6A12794E757CA78A, 0x00636E7C90AEFCCA},
     {0x90FDE64F18756308, 0xF8523408AE94B20D, 0x72210C61CF9F900A, 0x60C52B5A11A33276, 0x0EE66C4A82F6A731, 0x79B0BCEE1C4CDA78, 0xE0B16FC3FBA2166F, 0x014915E887177D09}},
    {{0x36D4DC7332E5B371, 0xC8093D457F7A631A, 0x5DB1B2034231F44B, 0xBC9BC856E8904F26, 0x6187B7E622BC2185, 0xF29D167A716A615E, 0x560A9338F538B46E, 0x008285CA01E612D4},
     {0xEF9CE312CA246A5B, 0xAB603CA21D4D4EA8, 0xC42B0D7E88BC630C, 0xF2A1E801F63EA496, 0xF6EB72E298DEF9EE, 0x5DCC7ED9B9D00D29, 0x8B9A274D2724999C, 0x017ABF685ACD835B}},
    {{0xF82DB2EDA77F028A, 0xA0B1FE3D2A190347, 0xAA52C70D2271D8D1, 0x5C2D69D6718E65D9, 0x60361BC498DF86C7, 0xF15A78DA16239358, 0x7E71872E93E2BDFD, 0x00E0CA9CAED594C1},
     {0x867D6D6EF66C887F, 0xB13BABC6178DB2A6, 0x1738B079A5BBE2B2, 0x3A9CD3358D593633, 0x083481ED9D6ECFEC, 0x5E7EE2DB93BC38AD, 0x90142DDC81E96226, 0x00EEFE4414788855}},
    {{0x3D0DF3BBA4744C92, 0xCDA70EF744DC5442, 0x15EDE5E5F2029614, 0x5DE88C7E9BF9618C, 0x41409BEE172E2A7E, 0x21234B85EF99A855, 0x57D4B80B9E9DC6BF, 0x00F1784B39DEA5BC},
     {0x84DD8D6957D0D49C, 0xAC8CA2F6A76A47BB, 0x11DBC66626E83704, 0x557360CE674147C0, 0x6733DB40F5488209, 0x153241627DC5627B, 0x3FC65D093FD65E2D, 0x0088C1DDC40561B1}},
    {{0x70837DB6E6B3A8BF, 0x892C4E3AB44F3A23, 0xF55BBB1A27151508, 0x57CFD8AB94ABD532, 0x2727C1B9E01C16CE, 0x0EF412F22F8A1550, 0xF9A1A4E4C90E2FBD, 0x0036D3B2C89FE2D5},
     {0xE10A41F013D57CBB, 0x4AF093241DA01790, 0x51C37561E08B7E7E, 0x513CFA8F62A89FE9, 0xCD9F0587A9C00813, 0x7192A0C83A51FA33, 0x764693A93B4F0924, 0x00571614255634B9}},
    {{0x1D16E6E000293CF0, 0x9F6D16F0E01BF2E8, 0x2538040B7A31CBA9, 0x5ECADBA94A37E39A, 0xBC47C74FD2047A82, 0xF37F9CFDC833849A, 0xF977EB94C9D067AE, 0x0107BAA101D19780},
     {0x69014EF349E1A56B, 0x1AC544F932E8295F, 0x77D09327740C06B4, 0x717BC9A4F75B9179, 0x8F7EC12983BBCD18, 0x2B246E75CE16D85D, 0x5B5BC2CD1947D558, 0x00B2E5544455AEA8}},
    {{0x857A30514395D836, 0x2C36EE9C7D66E8C7, 0x59081F50B68BC9B7, 0x53A653331BCDD8C8, 0xD5A1D5C73CB6FD46, 0xB3A0A34C7D1534D4, 0x13B52541D2FF6001, 0x01AFB4AD9CCBB922},
     {0xE00EB4EC26258B07, 0x2A549CD738AEE0EC, 0xF4763D9601E03988, 0xCB93FF6293F17595, 0x661D1D4B326C558F, 0x9948206ED0E14FC2, 0x881DBD71F3F042FF, 0x00205CB0E3F9A007}},
    {{0x17A639ABB543260A, 0x466B950EA32D242A, 0x1208BEE5EB968EE5, 0x1100E1944D861A90, 0xCD1C4C7EDD4A3756, 0x6887EE7D3186E7AF, 0x7288D3A522540C90, 0x008BDC6DE07EC8E6},
     {0x32E711F4108BD2B2, 0xF61CC8E4671BAED7, 0xA8516DB2F814C132, 0xCEA27FF21720A484, 0xB23B077AA25A9E71, 0xE53D7FB0441F5CC0, 0x4B9C7727CA4B4A9E, 0x011CDAA77B2403E4}},
    {{0xDB65B8E3B482CBE9, 0x6C9C9C698BC9A8AF, 0x73BF3C375A047947, 0x7F03CBAB7FEEE273, 0x9ABC6D69CE25D306, 0x431F2E910AD10853, 0x6BB5B7F5163EFC4A, 0x007E62B145A9B446},
     {0x54B6A71525A1897E, 0x52A8411DDDF47A5A, 0x295DF247093171DF, 0x1D271190C889D73F, 0xC955B002EE2F8C97, 0xF1C960A15514BB3F, 0x022EA8B55588465B, 0x01A5DF70CEDD6EC3}},
    {{0x5D97E84F9147A6D4, 0x26DF7BBC9C1AB182, 0xDC01C6970AB9A0C7, 0x7AF34A026BA266F5, 0xE9227C3F3318ECD4, 0x655AD013ABD17DE9, 0xA33A737686B57F91, 0x0074303BF0D54AC7},
     {0xBC1C86624806659B, 0x5E93399A6ADEC7CC, 0x76876F5B9F9B6AE8, 0x6ED31C1BFE10F2FD, 0x1CC27B731C2E9357, 0x7E9F116EE5192ABF, 0x9ACF6918E1874D22, 0x01801E640504CC76}},
    {{0x979256269B26FC58, 0x191CD184F177C18C, 0xC09715D3023BCBB0, 0xE770DE3644C91BF2, 0x0D388C4D335755BE, 0x351C94948B2E1BAB, 0x7792C8E76FAA3042, 0x00B45F87081C67F5},
     {0xDADC66D9A679C668, 0x4086F15E014AADCC, 0xC4BFB4A1DCB47B55, 0x79ACE16FCBE4AE72, 0x8A0865B0707E9EB0, 0xF51DE468EF80A5C2, 0x62A69FA320CCEDFE, 0x00DB0544EC1A4B2D}},
};

const uint64_t Z_NQR_TABLE[20][2][NWORDS_FIELD] = {
    {{0x0BD33F844011A886, 0x32040B6BD1D4A9C8, 0xE6D2F7C57CCC99C2, 0xC475F25D2FCB9D4D, 0x5B38621D6C71DEDB, 0x825B77BC8CF318DA, 0x8C1BFB53FDC60A95, 0x01645280E9E85449},
     {0xD96592793A6F5012, 0xD24E180DE58F0667, 0x2FF247F1D11558F2, 0x54488C648A7C1652, 0x00C0F1A178B97922, 0x75E2F682840603CB, 0xECCC4DCC0B0B7C14, 0x0109317FB221EF1D}},
    {{0xA7F6DDAD377E949F, 0xC9A6DA52DDB41A08, 0x2894672F2C697625, 0x8F4B442F7A2469B5, 0x7D6BAEC3AE6226DB, 0x597E5511964D0BFC, 0x7E834A5FD126D8D3, 0x0121DAFEF351CDA9},
     {0x153F526160FA3A68, 0xAC77E78E37787015, 0x817507796F95C726, 0x7190151DC1CC5E90, 0xB23567666204C6F0, 0xB0AE44BEDA7DEBC5, 0xBC13F503E4E448AA, 0x001FA2F31FB59B2F}},
    {{0xF19C43B88CDD60C9, 0x7202A06E7510ABD8, 0x83C2330AF6B3C229, 0xF2031F87D0303B42, 0x1330CBD252F82178, 0xA6E768DB07BB6A60, 0x0AB787C063435469, 0x00A51B46438119D6},
     {0x5A435C6708EFF376, 0xC70E23FC71F0A1DA, 0x64E002AEE87AFDD7, 0xF58E3EFCD6EB7ADA, 0x01DE298614E968CF, 0xF31F1589A59A91A7, 0x4A6CFCC45119C0B5, 0x00C518A6E3662477}},
    {{0x074150FA6F4D8E30, 0x5F386CA63D4885B7, 0xF2B7A89FAB264FB7, 0x76ECB7D18D69BA95, 0x7FD6771CD5DB5957, 0x2EFC9752C9B8133F, 0x4D0FD4B9A26B8498, 0x019989D71AD21C5E},
     {0x538902B24603B627, 0xF8B333DC6E4568E6, 0x8E309F4696D6A873, 0x5220B88563B8A110, 0x2CF33EDD0ABA96F3, 0x15C6D86286707189, 0x73AF589A46A9928E, 0x004A2BB2B0697D2A}},
    {{0xCB38EEB65FE708E0, 0xB5EE5EF196EE7A03, 0xB46DB052F99EF71D, 0xDA593D4F490FEC6E, 0xB625F30E389E96CC, 0x1E582359F90B2DEF, 0x2C95BE102ADC9D26, 0x00D71C545AC12808},
     {0x46E27AC36DA06760, 0x843DA216146C4428, 0xCEC64DBF7DBC67BB, 0x7CFC4718819ABC45, 0x526F4C613F4F82DD, 0x5A9254BFF7F15124, 0x5DEF2FAC260C74F8, 0x008388C5A65B8DA1}},
    {{0x9DE76EA6BA5E62D4, 0xF22B4F3C2883F1DA, 0xBBF58CB4A08BB9D0, 0x290B84D40F41CFE3, 0xCB5A1750CE5CB69E, 0x3771F3CE49B508ED, 0x4796E9A1DD5A9924, 0x01785EA89B9527E3},
     {0xFFAA42871309B614, 0xB73B1BB2592A8617, 0x8543D49CD8907A6D, 0x1F83CE8D8222F70D, 0x0ACC04E571A123A4, 0xFED2EB8CB33F6B86, 0xB12ED6AD8844B2A1, 0x0049EE9EE0B97EF6}},
    {{0x0F37AD561A33F162, 0x1A522CAF472366B2, 0xAAB4753B3EFCE21D, 0xFA1EAF5CF9B3608C, 0xCED46E93E38728E8, 0xD32E22A0F2C86108, 0x70D73D997C685D5E, 0x0089A64FF99DADD4},
     {0x2FC526C2E670C4E3, 0xB0AE0FDA50B72637, 0xB4BEDECAD2B25538, 0xB3F0442667686C04, 0xEC432C71DAB5A42B, 0xCD4BC043B1468DB3, 0x481A8F310A122833, 0x0168B26104997B37}},
    {{0x554A437E7674534C, 0x6F9410B93EC48793, 0x776BC9530DEAD1F8, 0x2C2D2B90C70B9C4E, 0xB2945898B024ADFE, 0xD6E9FF599A265313, 0xE3CB3DF015D4FBA0, 0x00888FF1C012D64B},
     {0xD923231A1BB2E9BD, 0x3142AB2E146C397A, 0xDA7A4DA00D6139D1, 0xB5BCD64A29C8B071, 0xB862A8E5AE8B77AD, 0x40ADE44F9D6A4609, 0x95884903F8C417F1, 0x0016E2EC4F5653CD}},
    {{0x94591A74FF7AB8A5, 0xEC4D28A6BE05A3FE, 0x1C73E358D7DAF97D, 0x3B76948F1F9D133B, 0xED4986B4AE832365, 0xE4F07777BBAA3581, 0xBCF73C2AC7EF08B9, 0x00E777FE57735CEB},
     {0x75CD80954C2E105B, 0x8E9E0188C5BDA4C2, 0x3A946F46BB53E84B, 0xA49630880177C62B, 0x955CF5CB828DFAF6, 0xCD04A34A37C872A9, 0x0329C98014AB5DAC, 0x008C5835FC4D8EEF}},
    {{0x269ECBEF0FB0D63D, 0x778D05252B67C323, 0x053C0022B4B7824B, 0x6088EFC818DECE2E, 0xCF5EC467E9704091, 0x492AAB926EF8DC83, 0x89D5A1D0902184C2, 0x00B38167F82E193B},
     {0x43C6AAC709CA46AA, 0xFDC9EAAA149A4AFF, 0x2633EF98231EDDAA, 0x05DEE463705DE278, 0xFFFD0B7BE8035CE6, 0x7F1624CC31AAC141, 0x0E4AA845D35DFA5C, 0x0152144AAB154ADB}},
    {{0x26C913E56AD80492, 0x410721AC3E625C2C, 0xFEB7756EFA5C2D5F, 0x14FCA9F9850FFD73, 0x3B322D42050BA9A6, 0xA9A0B8D05E8AAE21, 0x221FD18C91E782FD, 0x00302E592928FEF4},
     {0x61A1FC3C5E25F770, 0x4AC8842196B76090, 0xD782ABF12509379F, 0x991DC6E4E7E9DFE1, 0xEE506AB36B78BFB9, 0xC402491E5EDC5F79, 0xEB2F8FF07AC3CA8A, 0x00D68FF601C52E06}},
    {{0xFE2EBA6F9FADBBEB, 0x7761A161E7AA3B71, 0x06FB2929157C4F35, 0x2C7154C628B438EF, 0xCD97EED5DD3ED299, 0xB0591B9973D12A10, 0x5F00A229491E314F, 0x0188EC27A41D4B94},
     {0x5A289D6CC29E122D, 0x64B8F4FEC96D5550, 0x96355B6ECAD11F06, 0xCA84960CD6E30F1E, 0xBDE9DD44137F6825, 0x7D84739F4F9671FB, 0x3E5F960687372D9C, 0x0172C6AA93C01C05}},
    {{0x9AAE5100042A023E, 0x2C872AD473D41C67, 0x212BF2DC4C8FA557, 0xA16C00B50E0E0DD5, 0x12296CF3FF40FE14, 0x0CDA5DDFB014FCEB, 0x11B010BB440E58EF, 0x0142549D0D68162E},
     {0xF2D054929A7E4724, 0x89FD0670955919A6, 0x953E6A3B77681EC4, 0x82A09CC605C9F77F, 0x4143D7BC2BD376D2, 0xF90F194D1D4E2555, 0x023D751951C2F9A4, 0x007D8398842833A6}},
    {{0x6DC4DBA6299A7FF3, 0xB436DC90B9CD1429, 0xF5DF87D2C2D6562D, 0xED843458A58F7D65, 0x5F98BDF4DF3C5F7C, 0xC0E68FE8C0E78730, 0xDFB69A13731ECD87, 0x018EA9E102325E10},
     {0x1D0B5DCBB6CE9A9B, 0x6275D972360EF920, 0x1077EB6FFBC2D5E2, 0x9E68945F08701FEC, 0x2914105744534396, 0x1ADE76446417E05C, 0x7A4396A32F05DDCF, 0x00164CD2C01B9017}},
    {{0x3E9BD38A4196C113, 0x8469B0FEC63739CD, 0x68A5973243882E90, 0x1E3CCDCAF321A404, 0x9B62F5776A179260, 0x7056614F6DB5FDD8, 0xC38F70B71514E610, 0x00FA11F3F4EF89AD},
     {0xAD459ADE506D5494, 0xBA0DD06484ED9170, 0x5DA01177F3A0D1BC, 0x2860CD2DD1FDF1A0, 0x60D4BF59DFF8FD51, 0x2F028FDCDC7A8D86, 0xF0AA825FF80542D4, 0x00F5193A31255F66}},
    {{0x10967DAF104CB64A, 0x2D30666B3A66225B, 0x1FA669C19100E341, 0xAB21A7A83CC09E83, 0x608B11A55D020DA7, 0x5B21C0F46B5C0185, 0xC47B94494D0D3C18, 0x0077581C411DFA84},
     {0xA7EE2FEFD5224BC9, 0x25C5752EA6E91F96, 0xB5322418C6F08D6F, 0x362447C78FD1152A, 0x9C18254AB9EE981D, 0x272AE336F2707097, 0x00569607537AA7D0, 0x00E0AAB2D25ADC26}},
    {{0x98E2AC0376935413, 0xC72F51DEF4FD32CF, 0xE6E4370B9EB8B9AF, 0x30E89E2C8D352B29, 0x6C03F9CBBB2CF8B6, 0x779109F74D8A9422, 0x9D9292488BA4E016, 0x017BCC435B9EDBF1},
     {0x3B6FDBE617FE57E7, 0xB1326C93C077D31D, 0xCA909595225FFD6F, 0x9D932F915D4921D0, 0xFD9D65E2ECD9D360, 0x1B61493CE660EF48, 0x1A5DCC68CF60B6ED, 0x0193FCAC5DF00C69}},
    {{0x37DEE41AAA2F5BF8, 0xB912CF7A0C6A07BA, 0xB9989F1054CCF924, 0xE002FB52883517B1, 0x6C05655AB02B17C9, 0x37AFF27CDA26DD53, 0x6213A1D3553D3E5E, 0x007DADB0D63147D9},
     {0x2206689991F0FE9A, 0x4172826B1E5DAF30, 0xE0E8FFB4E6DED518, 0xA399075221B2200C, 0x6D1275097E83D6FF, 0x55152659BDA3E385, 0x2E83A36BE09A006C, 0x01230ED4BE5E283D}},
    {{0xBD66588F21BDA652, 0x9718DA3C8D197C33, 0xAA7C8BE9CCE1E870, 0x88373AA4055EC395, 0x8B4520657F9B2339, 0xD51E57BAC2DD5284, 0x1629EC3DC760D0FE, 0x01AFB228615A8F30},
     {0x3FD765040B767A4A, 0x11425FDE675FF266, 0x1C12017BBFEC6557, 0x421A07D78740FACF, 0x9D6B9C2EBD3B6880, 0xEB82FA1382AC112C, 0x0B3F24195A8D5B24, 0x00E6890BEB545FA8}},
    {{0xFE57B6C4FE44B30E, 0x8F0106D3C269F797, 0xCE4179F71273280A, 0xBE3016E70E5CF8F5, 0xB5EFDC2F8BD32DC9, 0x237E5219BADE8CB7, 0x82B746324716B6A8, 0x0138D178AEFB3C5A},
     {0x473EACBBC3B56D2F, 0x3BB269095A914A62, 0xA6B276A6FB61E18C, 0x91D93E303EF667A9, 0x40887540A26FBC42, 0xBA984BFAF31DE383, 0x91644F268C1ABF79, 0x00A84FB321E10361}},
};
