function evalDYmodDu(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	c144, c432, c395, c513, c466, c392, c91, c165, c263, c41, c416, c257, c511, c282, c504, c58, c2, c163, c128, c227, c453, c310, c320, c487, c12, c307, c390, c166, c96, c371, c175, c81, c529, c122, c495, c44, c530, c182, c80, c223, c302, c37, c289, c526, c255, c102, c425, c281, c501, c348, c158, c189, c366, c305, c116, c190, c42, c457, c410, c493, c76, c463, c444, c38, c17, c420, c210, c159, c385, c439, c370, c336, c11, c226, c314, c450, c127, c219, c98, c247, c64, c101, c315, c330, c482, c267, c142, c343, c123, c69, c60, c169, c253, c404, c552, c104, c114, c203, c94, c520, c36, c434, c364, c243, c555, c99, c424, c360, c156, c157, c313, c181, c440, c18, c264, c131, c489, c218, c300, c108, c137, c507, c95, c232, c465, c321, c229, c86, c377, c222, c230, c347, c344, c149, c522, c106, c280, c214, c550, c396, c409, c211, c57, c215, c7, c306, c151, c82, c148, c557, c361, c49, c473, c388, c92, c231, c192, c496, c153, c309, c237, c356, c183, c521, c475, c115, c193, c294, c452, c77, c207, c217, c422, c202, c10, c71, c462, c252, c63, c23, c55, c296, c30, c20, c54, c111, c331, c349, c113, c411, c204, c259, c239, c171, c240, c368, c248, c418, c194, c464, c312, c242, c297, c431, c415, c47, c133, c168, c34, c249, c25, c359, c503, c398, c502, c24, c523, c532, c3, c428, c209, c167, c185, c427, c332, c436, c379, c403, c235, c394, c33, c225, c383, c83, c105, c87, c442, c245, c455, c191, c468, c206, c363, c152, c279, c59, c201, c556, c558, c476, c346, c471, c399, c488, c393, c419, c265, c79, c467, c262, c484, c456, c75, c486, c129, c525, c384, c402, c246, c221, c417, c491, c50, c208, c236, c376, c266, c6, c407, c510, c103, c244, c479, c132, c375, c492, c519, c85, c499, c485, c299, c446, c500, c241, c46, c173, c89, c517, c387, c372, c172, c353, c480, c380, c469, c483, c48, c19, c195, c381, c213, c160, c308, c391, c170, c317, c90, c516, c139, c136, c303, c9, c39, c373, c260, c508, c62, c126, c197, c21, c506, c426, c518, c180, c31, c84, c65, c408, c340, c497, c35, c401, c14, c429, c205, c135, c386, c220, c472, c150, c196, c358, c478, c27, c357, c254, c286, c460, c430, c45, c147, c284, c494, c291, c339, c354, c216, c28, c400, c319, c70, c72, c8, c93, c13, c112, c198, c51, c67, c74, c16, c138, c22, c323, c338, c367, c382, c268, c141, c124, c448, c447, c374, c509, c52, c301, c551, c228, c88, c362, c125, c490, c421, c61, c365, c162, c68, c474, c130, c66, c412, c161, c324, c164, c15, c528, c389, c224, c449, c435, c337, c154, c143, c481, c311, c26, c184, c355, c531, c117, c341, c53, c512, c4, c405, c73, c145, c199, c238, c350, c233, c5, c454, c438, c287, c97, c258, c304, c250, c283, c335, c146, c269, c288, c322, c461, c351, c176, c325, c251, c56, c316, c451, c369, c78, c212, c174, c477, c200, c121, c437, c261, c445, c100, c470, c109, c553, c397, c290, c234, c554, c293, c140, c524, c155, c134, c527, c378, c505, c40, c318, c32, c256, c43, c498, c107, c345, c29, c110, c433= u 
	k71, kd25, kd120, c515, k120, kd47, kd94, k73, kd36, kd123h, kd119, kd49, kd4, k45, kd106b, kd10, kd71, k105, kd101, kd17, k106b, kd111, StepMini_bool3, k2, kd109, k111, kd23, k8, Step33_bool3, kd75, Step10Mio_bool2, kd60b, kd63, kd69, k94b, kd42, Step33_bool1, kd68b, k114, kd40, k74, kd114, k67, k69, k44, k94, kd33, kd5, Step10Mio_bool3, kd7, kd105, kd18, k5, k96, k76, k40, k102, kd35, k61, Step1_bool1, k117, kd56, Step5_bool1, kd104, kd115, k5b, kd50, k113, k20, c1, k112, kd5b, k35, k104, k75, k19, c514, kd97, k22, kd29, k65, kd2b, kd37, Step5_bool2, k50, kd76, k48, k25, kd8, k29, kd32, kd110, kd52, kd22, k107, k56, k60, k33, kd15, StepMini_bool1, kd67, Step1_bool3, k72, k58, kd108, k118, kd99, kd58, k101, k60b, k7, k34, k18, kd24, StepLate_bool1, kd34, kd112, k95, kd2, k53, kd103, kd117, k42, kd43, k70, k57, kd21, k6, Step10Mio_bool1, kd95, kd41, kd66, kd72, kd48, StepMini_bool2, Step0, k52, kd98, kd57, k49, k6b, k108, k41, kd60, kd64, k47, kd6b, k2b, kd55, kd107, k4, cyt, k115, k68, k109, kd61, k10b, k36, kd70, kd44, kd45, k16, k43, k66, k28, kd116, kd74, kd28, k64, Step33_bool2, kd73, k123, kd106, kd65, StepLate_bool2, k120b, c285, k103, kd118, kd97c, k62b, kd100, k37, kd68, StepLate_bool3, k110, k32, Step5_bool3, k17, k21, k4b, kd102, kd113, Step1_bool2, k15, k1d, k123h, kd53, kd22b, kd20, k23, k60c, k106, kd8b, kd123, kd19, k8b, kd96, k55 = p 
	if observableId == "AKT_PP" 
		observableParameter1_AKT_PP = getObsOrSdParam(obsPar, mapObsParam)
		out[261] = 1 / observableParameter1_AKT_PP
		out[340] = 1 / observableParameter1_AKT_PP
		out[349] = 1 / observableParameter1_AKT_PP
		out[495] = 1 / observableParameter1_AKT_PP
		return nothing
	end

	if observableId == "ERB_B1_P_tot" 
		observableParameter1_ERB_B1_P_tot = getObsOrSdParam(obsPar, mapObsParam)
		out[2] = 2 / observableParameter1_ERB_B1_P_tot
		out[5] = 1 / observableParameter1_ERB_B1_P_tot
		out[7] = 2 / observableParameter1_ERB_B1_P_tot
		out[8] = 1 / observableParameter1_ERB_B1_P_tot
		out[9] = 1 / observableParameter1_ERB_B1_P_tot
		out[11] = 1 / observableParameter1_ERB_B1_P_tot
		out[12] = 1 / observableParameter1_ERB_B1_P_tot
		out[14] = 1 / observableParameter1_ERB_B1_P_tot
		out[18] = 1 / observableParameter1_ERB_B1_P_tot
		out[20] = 1 / observableParameter1_ERB_B1_P_tot
		out[28] = 1 / observableParameter1_ERB_B1_P_tot
		out[29] = 2 / observableParameter1_ERB_B1_P_tot
		out[31] = 1 / observableParameter1_ERB_B1_P_tot
		out[38] = 1 / observableParameter1_ERB_B1_P_tot
		out[40] = 1 / observableParameter1_ERB_B1_P_tot
		out[42] = 2 / observableParameter1_ERB_B1_P_tot
		out[45] = 1 / observableParameter1_ERB_B1_P_tot
		out[48] = 1 / observableParameter1_ERB_B1_P_tot
		out[52] = 1 / observableParameter1_ERB_B1_P_tot
		out[56] = 1 / observableParameter1_ERB_B1_P_tot
		out[59] = 1 / observableParameter1_ERB_B1_P_tot
		out[65] = 2 / observableParameter1_ERB_B1_P_tot
		out[66] = 2 / observableParameter1_ERB_B1_P_tot
		out[67] = 1 / observableParameter1_ERB_B1_P_tot
		out[74] = 1 / observableParameter1_ERB_B1_P_tot
		out[76] = 1 / observableParameter1_ERB_B1_P_tot
		out[78] = 1 / observableParameter1_ERB_B1_P_tot
		out[79] = 2 / observableParameter1_ERB_B1_P_tot
		out[80] = 1 / observableParameter1_ERB_B1_P_tot
		out[81] = 2 / observableParameter1_ERB_B1_P_tot
		out[86] = 1 / observableParameter1_ERB_B1_P_tot
		out[93] = 1 / observableParameter1_ERB_B1_P_tot
		out[96] = 2 / observableParameter1_ERB_B1_P_tot
		out[98] = 1 / observableParameter1_ERB_B1_P_tot
		out[99] = 2 / observableParameter1_ERB_B1_P_tot
		out[101] = 2 / observableParameter1_ERB_B1_P_tot
		out[102] = 1 / observableParameter1_ERB_B1_P_tot
		out[104] = 1 / observableParameter1_ERB_B1_P_tot
		out[106] = 2 / observableParameter1_ERB_B1_P_tot
		out[112] = 1 / observableParameter1_ERB_B1_P_tot
		out[113] = 1 / observableParameter1_ERB_B1_P_tot
		out[114] = 2 / observableParameter1_ERB_B1_P_tot
		out[115] = 2 / observableParameter1_ERB_B1_P_tot
		out[116] = 1 / observableParameter1_ERB_B1_P_tot
		out[117] = 2 / observableParameter1_ERB_B1_P_tot
		out[118] = 1 / observableParameter1_ERB_B1_P_tot
		out[123] = 2 / observableParameter1_ERB_B1_P_tot
		out[124] = 1 / observableParameter1_ERB_B1_P_tot
		out[125] = 1 / observableParameter1_ERB_B1_P_tot
		out[127] = 1 / observableParameter1_ERB_B1_P_tot
		out[130] = 1 / observableParameter1_ERB_B1_P_tot
		out[131] = 1 / observableParameter1_ERB_B1_P_tot
		out[134] = 1 / observableParameter1_ERB_B1_P_tot
		out[138] = 1 / observableParameter1_ERB_B1_P_tot
		out[141] = 1 / observableParameter1_ERB_B1_P_tot
		out[142] = 1 / observableParameter1_ERB_B1_P_tot
		out[144] = 1 / observableParameter1_ERB_B1_P_tot
		out[145] = 2 / observableParameter1_ERB_B1_P_tot
		out[147] = 1 / observableParameter1_ERB_B1_P_tot
		out[149] = 1 / observableParameter1_ERB_B1_P_tot
		out[155] = 2 / observableParameter1_ERB_B1_P_tot
		out[156] = 1 / observableParameter1_ERB_B1_P_tot
		out[157] = 1 / observableParameter1_ERB_B1_P_tot
		out[159] = 1 / observableParameter1_ERB_B1_P_tot
		out[161] = 1 / observableParameter1_ERB_B1_P_tot
		out[163] = 1 / observableParameter1_ERB_B1_P_tot
		out[167] = 1 / observableParameter1_ERB_B1_P_tot
		out[171] = 1 / observableParameter1_ERB_B1_P_tot
		out[172] = 1 / observableParameter1_ERB_B1_P_tot
		out[174] = 1 / observableParameter1_ERB_B1_P_tot
		out[178] = 1 / observableParameter1_ERB_B1_P_tot
		out[179] = 2 / observableParameter1_ERB_B1_P_tot
		out[180] = 2 / observableParameter1_ERB_B1_P_tot
		out[184] = 2 / observableParameter1_ERB_B1_P_tot
		out[190] = 1 / observableParameter1_ERB_B1_P_tot
		out[191] = 1 / observableParameter1_ERB_B1_P_tot
		out[192] = 1 / observableParameter1_ERB_B1_P_tot
		out[193] = 1 / observableParameter1_ERB_B1_P_tot
		out[194] = 1 / observableParameter1_ERB_B1_P_tot
		out[195] = 1 / observableParameter1_ERB_B1_P_tot
		out[197] = 1 / observableParameter1_ERB_B1_P_tot
		out[199] = 1 / observableParameter1_ERB_B1_P_tot
		out[200] = 1 / observableParameter1_ERB_B1_P_tot
		out[202] = 1 / observableParameter1_ERB_B1_P_tot
		out[204] = 2 / observableParameter1_ERB_B1_P_tot
		out[205] = 2 / observableParameter1_ERB_B1_P_tot
		out[209] = 2 / observableParameter1_ERB_B1_P_tot
		out[210] = 1 / observableParameter1_ERB_B1_P_tot
		out[211] = 2 / observableParameter1_ERB_B1_P_tot
		out[220] = 1 / observableParameter1_ERB_B1_P_tot
		out[221] = 1 / observableParameter1_ERB_B1_P_tot
		out[222] = 1 / observableParameter1_ERB_B1_P_tot
		out[223] = 1 / observableParameter1_ERB_B1_P_tot
		out[224] = 1 / observableParameter1_ERB_B1_P_tot
		out[229] = 1 / observableParameter1_ERB_B1_P_tot
		out[231] = 2 / observableParameter1_ERB_B1_P_tot
		out[232] = 1 / observableParameter1_ERB_B1_P_tot
		out[238] = 1 / observableParameter1_ERB_B1_P_tot
		out[240] = 1 / observableParameter1_ERB_B1_P_tot
		out[242] = 1 / observableParameter1_ERB_B1_P_tot
		out[244] = 1 / observableParameter1_ERB_B1_P_tot
		out[247] = 1 / observableParameter1_ERB_B1_P_tot
		out[256] = 2 / observableParameter1_ERB_B1_P_tot
		out[257] = 1 / observableParameter1_ERB_B1_P_tot
		out[260] = 1 / observableParameter1_ERB_B1_P_tot
		out[264] = 2 / observableParameter1_ERB_B1_P_tot
		out[269] = 1 / observableParameter1_ERB_B1_P_tot
		out[270] = 1 / observableParameter1_ERB_B1_P_tot
		out[274] = 1 / observableParameter1_ERB_B1_P_tot
		out[275] = 1 / observableParameter1_ERB_B1_P_tot
		out[277] = 1 / observableParameter1_ERB_B1_P_tot
		out[282] = 1 / observableParameter1_ERB_B1_P_tot
		out[284] = 1 / observableParameter1_ERB_B1_P_tot
		out[292] = 1 / observableParameter1_ERB_B1_P_tot
		out[294] = 1 / observableParameter1_ERB_B1_P_tot
		out[296] = 1 / observableParameter1_ERB_B1_P_tot
		out[297] = 2 / observableParameter1_ERB_B1_P_tot
		out[301] = 1 / observableParameter1_ERB_B1_P_tot
		out[306] = 2 / observableParameter1_ERB_B1_P_tot
		out[308] = 2 / observableParameter1_ERB_B1_P_tot
		out[309] = 1 / observableParameter1_ERB_B1_P_tot
		out[311] = 1 / observableParameter1_ERB_B1_P_tot
		out[317] = 2 / observableParameter1_ERB_B1_P_tot
		out[320] = 2 / observableParameter1_ERB_B1_P_tot
		out[325] = 1 / observableParameter1_ERB_B1_P_tot
		out[329] = 1 / observableParameter1_ERB_B1_P_tot
		out[330] = 2 / observableParameter1_ERB_B1_P_tot
		out[334] = 1 / observableParameter1_ERB_B1_P_tot
		out[337] = 2 / observableParameter1_ERB_B1_P_tot
		out[341] = 2 / observableParameter1_ERB_B1_P_tot
		out[344] = 1 / observableParameter1_ERB_B1_P_tot
		out[345] = 1 / observableParameter1_ERB_B1_P_tot
		out[348] = 1 / observableParameter1_ERB_B1_P_tot
		out[350] = 1 / observableParameter1_ERB_B1_P_tot
		out[351] = 1 / observableParameter1_ERB_B1_P_tot
		out[354] = 2 / observableParameter1_ERB_B1_P_tot
		out[356] = 1 / observableParameter1_ERB_B1_P_tot
		out[367] = 1 / observableParameter1_ERB_B1_P_tot
		out[373] = 2 / observableParameter1_ERB_B1_P_tot
		out[374] = 2 / observableParameter1_ERB_B1_P_tot
		out[377] = 1 / observableParameter1_ERB_B1_P_tot
		out[379] = 2 / observableParameter1_ERB_B1_P_tot
		out[391] = 2 / observableParameter1_ERB_B1_P_tot
		out[392] = 1 / observableParameter1_ERB_B1_P_tot
		out[398] = 1 / observableParameter1_ERB_B1_P_tot
		out[399] = 2 / observableParameter1_ERB_B1_P_tot
		out[406] = 1 / observableParameter1_ERB_B1_P_tot
		out[407] = 2 / observableParameter1_ERB_B1_P_tot
		out[409] = 1 / observableParameter1_ERB_B1_P_tot
		out[410] = 2 / observableParameter1_ERB_B1_P_tot
		out[411] = 1 / observableParameter1_ERB_B1_P_tot
		out[414] = 1 / observableParameter1_ERB_B1_P_tot
		out[415] = 2 / observableParameter1_ERB_B1_P_tot
		out[418] = 1 / observableParameter1_ERB_B1_P_tot
		out[419] = 1 / observableParameter1_ERB_B1_P_tot
		out[420] = 1 / observableParameter1_ERB_B1_P_tot
		out[427] = 1 / observableParameter1_ERB_B1_P_tot
		out[438] = 1 / observableParameter1_ERB_B1_P_tot
		out[439] = 1 / observableParameter1_ERB_B1_P_tot
		out[441] = 1 / observableParameter1_ERB_B1_P_tot
		out[442] = 2 / observableParameter1_ERB_B1_P_tot
		out[444] = 1 / observableParameter1_ERB_B1_P_tot
		out[446] = 2 / observableParameter1_ERB_B1_P_tot
		out[447] = 1 / observableParameter1_ERB_B1_P_tot
		out[449] = 1 / observableParameter1_ERB_B1_P_tot
		out[458] = 1 / observableParameter1_ERB_B1_P_tot
		out[460] = 1 / observableParameter1_ERB_B1_P_tot
		out[463] = 1 / observableParameter1_ERB_B1_P_tot
		out[466] = 1 / observableParameter1_ERB_B1_P_tot
		out[467] = 1 / observableParameter1_ERB_B1_P_tot
		out[469] = 1 / observableParameter1_ERB_B1_P_tot
		out[471] = 1 / observableParameter1_ERB_B1_P_tot
		out[472] = 1 / observableParameter1_ERB_B1_P_tot
		out[473] = 1 / observableParameter1_ERB_B1_P_tot
		out[474] = 2 / observableParameter1_ERB_B1_P_tot
		out[480] = 1 / observableParameter1_ERB_B1_P_tot
		out[492] = 2 / observableParameter1_ERB_B1_P_tot
		out[493] = 1 / observableParameter1_ERB_B1_P_tot
		out[498] = 2 / observableParameter1_ERB_B1_P_tot
		out[500] = 1 / observableParameter1_ERB_B1_P_tot
		return nothing
	end

	if observableId == "ERK_PP" 
		observableParameter1_ERK_PP = getObsOrSdParam(obsPar, mapObsParam)
		out[2] = 1 / observableParameter1_ERK_PP
		out[29] = 1 / observableParameter1_ERK_PP
		out[46] = 1 / observableParameter1_ERK_PP
		out[79] = 1 / observableParameter1_ERK_PP
		out[82] = 1 / observableParameter1_ERK_PP
		out[102] = 1 / observableParameter1_ERK_PP
		out[113] = 1 / observableParameter1_ERK_PP
		out[123] = 1 / observableParameter1_ERK_PP
		out[165] = 1 / observableParameter1_ERK_PP
		out[204] = 1 / observableParameter1_ERK_PP
		out[234] = 1 / observableParameter1_ERK_PP
		out[246] = 1 / observableParameter1_ERK_PP
		out[303] = 1 / observableParameter1_ERK_PP
		out[336] = 1 / observableParameter1_ERK_PP
		out[353] = 1 / observableParameter1_ERK_PP
		out[404] = 1 / observableParameter1_ERK_PP
		out[408] = 1 / observableParameter1_ERK_PP
		out[420] = 1 / observableParameter1_ERK_PP
		out[424] = 1 / observableParameter1_ERK_PP
		out[444] = 1 / observableParameter1_ERK_PP
		out[446] = 1 / observableParameter1_ERK_PP
		out[468] = 1 / observableParameter1_ERK_PP
		out[471] = 1 / observableParameter1_ERK_PP
		out[500] = 1 / observableParameter1_ERK_PP
		return nothing
	end

end

function evalDYmodDp(u, t, p, obsPar, nonDynParam, observableId, mapObsParam, out) 

	c144, c432, c395, c513, c466, c392, c91, c165, c263, c41, c416, c257, c511, c282, c504, c58, c2, c163, c128, c227, c453, c310, c320, c487, c12, c307, c390, c166, c96, c371, c175, c81, c529, c122, c495, c44, c530, c182, c80, c223, c302, c37, c289, c526, c255, c102, c425, c281, c501, c348, c158, c189, c366, c305, c116, c190, c42, c457, c410, c493, c76, c463, c444, c38, c17, c420, c210, c159, c385, c439, c370, c336, c11, c226, c314, c450, c127, c219, c98, c247, c64, c101, c315, c330, c482, c267, c142, c343, c123, c69, c60, c169, c253, c404, c552, c104, c114, c203, c94, c520, c36, c434, c364, c243, c555, c99, c424, c360, c156, c157, c313, c181, c440, c18, c264, c131, c489, c218, c300, c108, c137, c507, c95, c232, c465, c321, c229, c86, c377, c222, c230, c347, c344, c149, c522, c106, c280, c214, c550, c396, c409, c211, c57, c215, c7, c306, c151, c82, c148, c557, c361, c49, c473, c388, c92, c231, c192, c496, c153, c309, c237, c356, c183, c521, c475, c115, c193, c294, c452, c77, c207, c217, c422, c202, c10, c71, c462, c252, c63, c23, c55, c296, c30, c20, c54, c111, c331, c349, c113, c411, c204, c259, c239, c171, c240, c368, c248, c418, c194, c464, c312, c242, c297, c431, c415, c47, c133, c168, c34, c249, c25, c359, c503, c398, c502, c24, c523, c532, c3, c428, c209, c167, c185, c427, c332, c436, c379, c403, c235, c394, c33, c225, c383, c83, c105, c87, c442, c245, c455, c191, c468, c206, c363, c152, c279, c59, c201, c556, c558, c476, c346, c471, c399, c488, c393, c419, c265, c79, c467, c262, c484, c456, c75, c486, c129, c525, c384, c402, c246, c221, c417, c491, c50, c208, c236, c376, c266, c6, c407, c510, c103, c244, c479, c132, c375, c492, c519, c85, c499, c485, c299, c446, c500, c241, c46, c173, c89, c517, c387, c372, c172, c353, c480, c380, c469, c483, c48, c19, c195, c381, c213, c160, c308, c391, c170, c317, c90, c516, c139, c136, c303, c9, c39, c373, c260, c508, c62, c126, c197, c21, c506, c426, c518, c180, c31, c84, c65, c408, c340, c497, c35, c401, c14, c429, c205, c135, c386, c220, c472, c150, c196, c358, c478, c27, c357, c254, c286, c460, c430, c45, c147, c284, c494, c291, c339, c354, c216, c28, c400, c319, c70, c72, c8, c93, c13, c112, c198, c51, c67, c74, c16, c138, c22, c323, c338, c367, c382, c268, c141, c124, c448, c447, c374, c509, c52, c301, c551, c228, c88, c362, c125, c490, c421, c61, c365, c162, c68, c474, c130, c66, c412, c161, c324, c164, c15, c528, c389, c224, c449, c435, c337, c154, c143, c481, c311, c26, c184, c355, c531, c117, c341, c53, c512, c4, c405, c73, c145, c199, c238, c350, c233, c5, c454, c438, c287, c97, c258, c304, c250, c283, c335, c146, c269, c288, c322, c461, c351, c176, c325, c251, c56, c316, c451, c369, c78, c212, c174, c477, c200, c121, c437, c261, c445, c100, c470, c109, c553, c397, c290, c234, c554, c293, c140, c524, c155, c134, c527, c378, c505, c40, c318, c32, c256, c43, c498, c107, c345, c29, c110, c433= u 
	k71, kd25, kd120, c515, k120, kd47, kd94, k73, kd36, kd123h, kd119, kd49, kd4, k45, kd106b, kd10, kd71, k105, kd101, kd17, k106b, kd111, StepMini_bool3, k2, kd109, k111, kd23, k8, Step33_bool3, kd75, Step10Mio_bool2, kd60b, kd63, kd69, k94b, kd42, Step33_bool1, kd68b, k114, kd40, k74, kd114, k67, k69, k44, k94, kd33, kd5, Step10Mio_bool3, kd7, kd105, kd18, k5, k96, k76, k40, k102, kd35, k61, Step1_bool1, k117, kd56, Step5_bool1, kd104, kd115, k5b, kd50, k113, k20, c1, k112, kd5b, k35, k104, k75, k19, c514, kd97, k22, kd29, k65, kd2b, kd37, Step5_bool2, k50, kd76, k48, k25, kd8, k29, kd32, kd110, kd52, kd22, k107, k56, k60, k33, kd15, StepMini_bool1, kd67, Step1_bool3, k72, k58, kd108, k118, kd99, kd58, k101, k60b, k7, k34, k18, kd24, StepLate_bool1, kd34, kd112, k95, kd2, k53, kd103, kd117, k42, kd43, k70, k57, kd21, k6, Step10Mio_bool1, kd95, kd41, kd66, kd72, kd48, StepMini_bool2, Step0, k52, kd98, kd57, k49, k6b, k108, k41, kd60, kd64, k47, kd6b, k2b, kd55, kd107, k4, cyt, k115, k68, k109, kd61, k10b, k36, kd70, kd44, kd45, k16, k43, k66, k28, kd116, kd74, kd28, k64, Step33_bool2, kd73, k123, kd106, kd65, StepLate_bool2, k120b, c285, k103, kd118, kd97c, k62b, kd100, k37, kd68, StepLate_bool3, k110, k32, Step5_bool3, k17, k21, k4b, kd102, kd113, Step1_bool2, k15, k1d, k123h, kd53, kd22b, kd20, k23, k60c, k106, kd8b, kd123, kd19, k8b, kd96, k55 = p 
	if observableId == "AKT_PP" 
		return nothing
	end

	if observableId == "ERB_B1_P_tot" 
		return nothing
	end

	if observableId == "ERK_PP" 
		return nothing
	end

end

function evalDSdDu!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	c144, c432, c395, c513, c466, c392, c91, c165, c263, c41, c416, c257, c511, c282, c504, c58, c2, c163, c128, c227, c453, c310, c320, c487, c12, c307, c390, c166, c96, c371, c175, c81, c529, c122, c495, c44, c530, c182, c80, c223, c302, c37, c289, c526, c255, c102, c425, c281, c501, c348, c158, c189, c366, c305, c116, c190, c42, c457, c410, c493, c76, c463, c444, c38, c17, c420, c210, c159, c385, c439, c370, c336, c11, c226, c314, c450, c127, c219, c98, c247, c64, c101, c315, c330, c482, c267, c142, c343, c123, c69, c60, c169, c253, c404, c552, c104, c114, c203, c94, c520, c36, c434, c364, c243, c555, c99, c424, c360, c156, c157, c313, c181, c440, c18, c264, c131, c489, c218, c300, c108, c137, c507, c95, c232, c465, c321, c229, c86, c377, c222, c230, c347, c344, c149, c522, c106, c280, c214, c550, c396, c409, c211, c57, c215, c7, c306, c151, c82, c148, c557, c361, c49, c473, c388, c92, c231, c192, c496, c153, c309, c237, c356, c183, c521, c475, c115, c193, c294, c452, c77, c207, c217, c422, c202, c10, c71, c462, c252, c63, c23, c55, c296, c30, c20, c54, c111, c331, c349, c113, c411, c204, c259, c239, c171, c240, c368, c248, c418, c194, c464, c312, c242, c297, c431, c415, c47, c133, c168, c34, c249, c25, c359, c503, c398, c502, c24, c523, c532, c3, c428, c209, c167, c185, c427, c332, c436, c379, c403, c235, c394, c33, c225, c383, c83, c105, c87, c442, c245, c455, c191, c468, c206, c363, c152, c279, c59, c201, c556, c558, c476, c346, c471, c399, c488, c393, c419, c265, c79, c467, c262, c484, c456, c75, c486, c129, c525, c384, c402, c246, c221, c417, c491, c50, c208, c236, c376, c266, c6, c407, c510, c103, c244, c479, c132, c375, c492, c519, c85, c499, c485, c299, c446, c500, c241, c46, c173, c89, c517, c387, c372, c172, c353, c480, c380, c469, c483, c48, c19, c195, c381, c213, c160, c308, c391, c170, c317, c90, c516, c139, c136, c303, c9, c39, c373, c260, c508, c62, c126, c197, c21, c506, c426, c518, c180, c31, c84, c65, c408, c340, c497, c35, c401, c14, c429, c205, c135, c386, c220, c472, c150, c196, c358, c478, c27, c357, c254, c286, c460, c430, c45, c147, c284, c494, c291, c339, c354, c216, c28, c400, c319, c70, c72, c8, c93, c13, c112, c198, c51, c67, c74, c16, c138, c22, c323, c338, c367, c382, c268, c141, c124, c448, c447, c374, c509, c52, c301, c551, c228, c88, c362, c125, c490, c421, c61, c365, c162, c68, c474, c130, c66, c412, c161, c324, c164, c15, c528, c389, c224, c449, c435, c337, c154, c143, c481, c311, c26, c184, c355, c531, c117, c341, c53, c512, c4, c405, c73, c145, c199, c238, c350, c233, c5, c454, c438, c287, c97, c258, c304, c250, c283, c335, c146, c269, c288, c322, c461, c351, c176, c325, c251, c56, c316, c451, c369, c78, c212, c174, c477, c200, c121, c437, c261, c445, c100, c470, c109, c553, c397, c290, c234, c554, c293, c140, c524, c155, c134, c527, c378, c505, c40, c318, c32, c256, c43, c498, c107, c345, c29, c110, c433= u 
	k71, kd25, kd120, c515, k120, kd47, kd94, k73, kd36, kd123h, kd119, kd49, kd4, k45, kd106b, kd10, kd71, k105, kd101, kd17, k106b, kd111, StepMini_bool3, k2, kd109, k111, kd23, k8, Step33_bool3, kd75, Step10Mio_bool2, kd60b, kd63, kd69, k94b, kd42, Step33_bool1, kd68b, k114, kd40, k74, kd114, k67, k69, k44, k94, kd33, kd5, Step10Mio_bool3, kd7, kd105, kd18, k5, k96, k76, k40, k102, kd35, k61, Step1_bool1, k117, kd56, Step5_bool1, kd104, kd115, k5b, kd50, k113, k20, c1, k112, kd5b, k35, k104, k75, k19, c514, kd97, k22, kd29, k65, kd2b, kd37, Step5_bool2, k50, kd76, k48, k25, kd8, k29, kd32, kd110, kd52, kd22, k107, k56, k60, k33, kd15, StepMini_bool1, kd67, Step1_bool3, k72, k58, kd108, k118, kd99, kd58, k101, k60b, k7, k34, k18, kd24, StepLate_bool1, kd34, kd112, k95, kd2, k53, kd103, kd117, k42, kd43, k70, k57, kd21, k6, Step10Mio_bool1, kd95, kd41, kd66, kd72, kd48, StepMini_bool2, Step0, k52, kd98, kd57, k49, k6b, k108, k41, kd60, kd64, k47, kd6b, k2b, kd55, kd107, k4, cyt, k115, k68, k109, kd61, k10b, k36, kd70, kd44, kd45, k16, k43, k66, k28, kd116, kd74, kd28, k64, Step33_bool2, kd73, k123, kd106, kd65, StepLate_bool2, k120b, c285, k103, kd118, kd97c, k62b, kd100, k37, kd68, StepLate_bool3, k110, k32, Step5_bool3, k17, k21, k4b, kd102, kd113, Step1_bool2, k15, k1d, k123h, kd53, kd22b, kd20, k23, k60c, k106, kd8b, kd123, kd19, k8b, kd96, k55 = p 
	if observableId == "AKT_PP" 
		return nothing
	end

	if observableId == "ERB_B1_P_tot" 
		return nothing
	end

	if observableId == "ERK_PP" 
		return nothing
	end

end

function evalDSdDp!(u, t, sdPar, p, nonDynParam, paramData, observableId, mapSdParam, out) 

	c144, c432, c395, c513, c466, c392, c91, c165, c263, c41, c416, c257, c511, c282, c504, c58, c2, c163, c128, c227, c453, c310, c320, c487, c12, c307, c390, c166, c96, c371, c175, c81, c529, c122, c495, c44, c530, c182, c80, c223, c302, c37, c289, c526, c255, c102, c425, c281, c501, c348, c158, c189, c366, c305, c116, c190, c42, c457, c410, c493, c76, c463, c444, c38, c17, c420, c210, c159, c385, c439, c370, c336, c11, c226, c314, c450, c127, c219, c98, c247, c64, c101, c315, c330, c482, c267, c142, c343, c123, c69, c60, c169, c253, c404, c552, c104, c114, c203, c94, c520, c36, c434, c364, c243, c555, c99, c424, c360, c156, c157, c313, c181, c440, c18, c264, c131, c489, c218, c300, c108, c137, c507, c95, c232, c465, c321, c229, c86, c377, c222, c230, c347, c344, c149, c522, c106, c280, c214, c550, c396, c409, c211, c57, c215, c7, c306, c151, c82, c148, c557, c361, c49, c473, c388, c92, c231, c192, c496, c153, c309, c237, c356, c183, c521, c475, c115, c193, c294, c452, c77, c207, c217, c422, c202, c10, c71, c462, c252, c63, c23, c55, c296, c30, c20, c54, c111, c331, c349, c113, c411, c204, c259, c239, c171, c240, c368, c248, c418, c194, c464, c312, c242, c297, c431, c415, c47, c133, c168, c34, c249, c25, c359, c503, c398, c502, c24, c523, c532, c3, c428, c209, c167, c185, c427, c332, c436, c379, c403, c235, c394, c33, c225, c383, c83, c105, c87, c442, c245, c455, c191, c468, c206, c363, c152, c279, c59, c201, c556, c558, c476, c346, c471, c399, c488, c393, c419, c265, c79, c467, c262, c484, c456, c75, c486, c129, c525, c384, c402, c246, c221, c417, c491, c50, c208, c236, c376, c266, c6, c407, c510, c103, c244, c479, c132, c375, c492, c519, c85, c499, c485, c299, c446, c500, c241, c46, c173, c89, c517, c387, c372, c172, c353, c480, c380, c469, c483, c48, c19, c195, c381, c213, c160, c308, c391, c170, c317, c90, c516, c139, c136, c303, c9, c39, c373, c260, c508, c62, c126, c197, c21, c506, c426, c518, c180, c31, c84, c65, c408, c340, c497, c35, c401, c14, c429, c205, c135, c386, c220, c472, c150, c196, c358, c478, c27, c357, c254, c286, c460, c430, c45, c147, c284, c494, c291, c339, c354, c216, c28, c400, c319, c70, c72, c8, c93, c13, c112, c198, c51, c67, c74, c16, c138, c22, c323, c338, c367, c382, c268, c141, c124, c448, c447, c374, c509, c52, c301, c551, c228, c88, c362, c125, c490, c421, c61, c365, c162, c68, c474, c130, c66, c412, c161, c324, c164, c15, c528, c389, c224, c449, c435, c337, c154, c143, c481, c311, c26, c184, c355, c531, c117, c341, c53, c512, c4, c405, c73, c145, c199, c238, c350, c233, c5, c454, c438, c287, c97, c258, c304, c250, c283, c335, c146, c269, c288, c322, c461, c351, c176, c325, c251, c56, c316, c451, c369, c78, c212, c174, c477, c200, c121, c437, c261, c445, c100, c470, c109, c553, c397, c290, c234, c554, c293, c140, c524, c155, c134, c527, c378, c505, c40, c318, c32, c256, c43, c498, c107, c345, c29, c110, c433= u 
	k71, kd25, kd120, c515, k120, kd47, kd94, k73, kd36, kd123h, kd119, kd49, kd4, k45, kd106b, kd10, kd71, k105, kd101, kd17, k106b, kd111, StepMini_bool3, k2, kd109, k111, kd23, k8, Step33_bool3, kd75, Step10Mio_bool2, kd60b, kd63, kd69, k94b, kd42, Step33_bool1, kd68b, k114, kd40, k74, kd114, k67, k69, k44, k94, kd33, kd5, Step10Mio_bool3, kd7, kd105, kd18, k5, k96, k76, k40, k102, kd35, k61, Step1_bool1, k117, kd56, Step5_bool1, kd104, kd115, k5b, kd50, k113, k20, c1, k112, kd5b, k35, k104, k75, k19, c514, kd97, k22, kd29, k65, kd2b, kd37, Step5_bool2, k50, kd76, k48, k25, kd8, k29, kd32, kd110, kd52, kd22, k107, k56, k60, k33, kd15, StepMini_bool1, kd67, Step1_bool3, k72, k58, kd108, k118, kd99, kd58, k101, k60b, k7, k34, k18, kd24, StepLate_bool1, kd34, kd112, k95, kd2, k53, kd103, kd117, k42, kd43, k70, k57, kd21, k6, Step10Mio_bool1, kd95, kd41, kd66, kd72, kd48, StepMini_bool2, Step0, k52, kd98, kd57, k49, k6b, k108, k41, kd60, kd64, k47, kd6b, k2b, kd55, kd107, k4, cyt, k115, k68, k109, kd61, k10b, k36, kd70, kd44, kd45, k16, k43, k66, k28, kd116, kd74, kd28, k64, Step33_bool2, kd73, k123, kd106, kd65, StepLate_bool2, k120b, c285, k103, kd118, kd97c, k62b, kd100, k37, kd68, StepLate_bool3, k110, k32, Step5_bool3, k17, k21, k4b, kd102, kd113, Step1_bool2, k15, k1d, k123h, kd53, kd22b, kd20, k23, k60c, k106, kd8b, kd123, kd19, k8b, kd96, k55 = p 
	if observableId == "AKT_PP" 
		return nothing
	end

	if observableId == "ERB_B1_P_tot" 
		return nothing
	end

	if observableId == "ERK_PP" 
		return nothing
	end

end

