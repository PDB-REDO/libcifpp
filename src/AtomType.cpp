// Lib for working with structures as contained in mmCIF and PDB files

#include "cif++/AtomType.h"
#include "cif++/Cif++.h"

using namespace std;

namespace libcif
{

const float kNA = nan("1");

const AtomTypeInfo kKnownAtoms[] =
{
	{ Nn,	"Unknown",			"Nn",	0,		false, {	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },
	{ H,	"Hydro­gen",			"H",	1.008,	false, {	53,		25,		37,		32,		kNA,	kNA,	120 } },
	{ He,	"He­lium",			"He",	4.0026,	false, {	31,		kNA,	32,		46,		kNA,	kNA,	140 } },
	{ Li,	"Lith­ium",			"Li",	6.94,	true,  {	167,	145,	134,	133,	124,	kNA,	182 } },
	{ Be,	"Beryl­lium",		"Be",	9.0122,	true,  {	112,	105,	90,		102,	90,		85,		kNA } },
	{ B,	"Boron",			"B",	10.81,	true,  {	87,		85,		82,		85,		78,		73,		kNA } },
	{ C,	"Carbon",			"C",	12.011,	false, {	67,		70,		77,		75,		67,		60,		170 } },
	{ N,	"Nitro­gen",			"N",	14.007,	false, {	56,		65,		75,		71,		60,		54,		155 } },
	{ O,	"Oxy­gen",			"O",	15.999,	false, {	48,		60,		73,		63,		57,		53,		152 } },
	{ F,	"Fluor­ine",			"F",	18.998,	false, {	42,		50,		71,		64,		59,		53,		147 } },
	{ Ne,	"Neon",				"Ne",	20.180,	false, {	38,		kNA,	69,		67,		96,		kNA,	154 } },
	{ Na,	"So­dium",			"Na",	22.990,	true,  {	190,	180,	154,	155,	160,	kNA,	227 } },
	{ Mg,	"Magne­sium",		"Mg",	24.305,	true,  {	145,	150,	130,	139,	132,	127,	173 } },
	{ Al,	"Alumin­ium",		"Al",	26.982,	true,  {	118,	125,	118,	126,	113,	111,	kNA } },
	{ Si,	"Sili­con",			"Si",	28.085,	true,  {	111,	110,	111,	116,	107,	102,	210 } },
	{ P,	"Phos­phorus",		"P",	30.974,	false, {	98,		100,	106,	111,	102,	94,		180 } },
	{ S,	"Sulfur",			"S",	32.06,	false, {	88,		100,	102,	103,	94,		95,		180 } },
	{ Cl,	"Chlor­ine",			"Cl",	35.45,	false, {	79,		100,	99,		99,		95,		93,		175 } },
	{ Ar,	"Argon",			"Ar",	39.948,	false, {	71,		kNA,	97,		96,		107,	96,		188 } },
	{ K,	"Potas­sium",		"K",	39.098,	true,  {	243,	220,	196,	196,	193,	kNA,	275 } },
	{ Ca,	"Cal­cium",			"Ca",	40.078,	true,  {	194,	180,	174,	171,	147,	133,	kNA } },
	{ Sc,	"Scan­dium",			"Sc",	44.956,	true,  {	184,	160,	144,	148,	116,	114,	kNA } },
	{ Ti,	"Tita­nium",			"Ti",	47.867,	true,  {	176,	140,	136,	136,	117,	108,	kNA } },
	{ V,	"Vana­dium",			"V",	50.942,	true,  {	171,	135,	125,	134,	112,	106,	kNA } },
	{ Cr,	"Chrom­ium",			"Cr",	51.996,	true,  {	166,	140,	127,	122,	111,	103,	kNA } },
	{ Mn,	"Manga­nese",		"Mn",	54.938,	true,  {	161,	140,	139,	119,	105,	103,	kNA } },
	{ Fe,	"Iron",				"Fe",	55.845,	true,  {	156,	140,	125,	116,	109,	102,	kNA } },
	{ Co,	"Cobalt",			"Co",	58.933,	true,  {	152,	135,	126,	111,	103,	96,		kNA } },
	{ Ni,	"Nickel",			"Ni",	58.693,	true,  {	149,	135,	121,	110,	101,	101,	163 } },
	{ Cu,	"Copper",			"Cu",	63.546,	true,  {	145,	135,	138,	112,	115,	120,	140 } },
	{ Zn,	"Zinc",				"Zn",	65.38,	true,  {	142,	135,	131,	118,	120,	kNA,	139 } },
	{ Ga,	"Gallium",			"Ga",	69.723,	true,  {	136,	130,	126,	124,	117,	121,	187 } },
	{ Ge,	"Germa­nium",		"Ge",	72.630,	true,  {	125,	125,	122,	121,	111,	114,	kNA } },
	{ As,	"Arsenic",			"As",	74.922,	true,  {	114,	115,	119,	121,	114,	106,	185 } },
	{ Se,	"Sele­nium",			"Se",	78.971,	false, {	103,	115,	116,	116,	107,	107,	190 } },
	{ Br,	"Bromine",			"Br",	79.904,	false, {	94,		115,	114,	114,	109,	110,	185 } },
	{ Kr,	"Kryp­ton",			"Kr",	83.798,	false, {	88,		kNA,	110,	117,	121,	108,	202 } },
	{ Rb,	"Rubid­ium",			"Rb",	85.468,	true,  {	265,	235,	211,	210,	202,	kNA,	kNA } },
	{ Sr,	"Stront­ium",		"Sr",	87.62,	true,  {	219,	200,	192,	185,	157,	139,	kNA } },
	{ Y,	"Yttrium",			"Y",	88.906,	true,  {	212,	180,	162,	163,	130,	124,	kNA } },
	{ Zr,	"Zirco­nium",		"Zr",	91.224,	true,  {	206,	155,	148,	154,	127,	121,	kNA } },
	{ Nb,	"Nio­bium",			"Nb",	92.906,	true,  {	198,	145,	137,	147,	125,	116,	kNA } },
	{ Mo,	"Molyb­denum",		"Mo",	95.95,	true,  {	190,	145,	145,	138,	121,	113,	kNA } },
	{ Tc,	"Tech­netium",		"Tc",	98,		true,  {	183,	135,	156,	128,	120,	110,	kNA } },
	{ Ru,	"Ruthe­nium",		"Ru",	101.07,	true,  {	178,	130,	126,	125,	114,	103,	kNA } },
	{ Rh,	"Rho­dium",			"Rh",	102.91,	true,  {	173,	135,	135,	125,	110,	106,	kNA } },
	{ Pd,	"Pallad­ium",		"Pd",	106.42,	true,  {	169,	140,	131,	120,	117,	112,	163 } },
	{ Ag,	"Silver",			"Ag",	107.87,	true,  {	165,	160,	153,	128,	139,	137,	172 } },
	{ Cd,	"Cad­mium",			"Cd",	112.41,	true,  {	161,	155,	148,	136,	144,	kNA,	158 } },
	{ In,	"Indium",			"In",	114.82,	true,  {	156,	155,	144,	142,	136,	146,	193 } },
	{ Sn,	"Tin",				"Sn",	118.71,	true,  {	145,	145,	141,	140,	130,	132,	217 } },
	{ Sb,	"Anti­mony",			"Sb",	121.76,	false, {	133,	145,	138,	140,	133,	127,	kNA } },
	{ Te,	"Tellurium",		"Te",	127.60,	false, {	123,	140,	135,	136,	128,	121,	206 } },
	{ I,	"Iodine",			"I",	126.90,	false, {	115,	140,	133,	133,	129,	125,	198 } },
	{ Xe,	"Xenon",			"Xe",	131.29,	false, {	108,	kNA,	130,	131,	135,	122,	216 } },
	{ Cs,	"Cae­sium",			"Cs",	132.91,	true,  {	298,	260,	225,	232,	209,	kNA,	kNA } },
	{ Ba,	"Ba­rium",			"Ba",	137.33,	true,  {	253,	215,	198,	196,	161,	149,	kNA } },
	{ La,	"Lan­thanum",		"La",	138.91,	true,  {	kNA,	195,	169,	180,	139,	139,	kNA } },
	{ Hf,	"Haf­nium",			"Hf",	178.49,	true,  {	208,	155,	150,	152,	128,	122,	kNA } },
	{ Ta,	"Tanta­lum",			"Ta",	180.95,	true,  {	200,	145,	138,	146,	126,	119,	kNA } },
	{ W,	"Tung­sten",			"W",	183.84,	true,  {	193,	135,	146,	137,	120,	115,	kNA } },
	{ Re,	"Rhe­nium",			"Re",	186.21,	true,  {	188,	135,	159,	131,	119,	110,	kNA } },
	{ Os,	"Os­mium",			"Os",	190.23,	true,  {	185,	130,	128,	129,	116,	109,	kNA } },
	{ Ir,	"Iridium",			"Ir",	192.22,	true,  {	180,	135,	137,	122,	115,	107,	kNA } },
	{ Pt,	"Plat­inum",			"Pt",	195.08,	true,  {	177,	135,	128,	123,	112,	110,	175 } },
	{ Au,	"Gold",				"Au",	196.97,	true,  {	174,	135,	144,	124,	121,	123,	166 } },
	{ Hg,	"Mer­cury",			"Hg",	200.59,	true,  {	171,	150,	149,	133,	142,	kNA,	155 } },
	{ Tl,	"Thallium",			"Tl",	204.38,	true,  {	156,	190,	148,	144,	142,	150,	196 } },
	{ Pb,	"Lead",				"Pb",	207.2,	true,  {	154,	180,	147,	144,	135,	137,	202 } },
	{ Bi,	"Bis­muth",			"Bi",	208.98,	true,  {	143,	160,	146,	151,	141,	135,	kNA } },
	{ Po,	"Polo­nium",			"Po",	209,	true,  {	135,	190,	kNA,	145,	135,	129,	kNA } },
	{ At,	"Asta­tine",			"At",	210,	false, {	127,	kNA,	kNA,	147,	138,	138,	kNA } },
	{ Rn,	"Radon",			"Rn",	222,	false, {	120,	kNA,	145,	142,	145,	133,	kNA } },
	{ Fr,	"Fran­cium",			"Fr",	223,	true,  {	kNA,	kNA,	kNA,	223,	218,	kNA,	kNA } },
	{ Ra,	"Ra­dium",			"Ra",	226,	true,  {	kNA,	215,	kNA,	201,	173,	159,	kNA } },
	{ Ac,	"Actin­ium",			"Ac",	227,	true,  {	kNA,	195,	kNA,	186,	153,	140,	kNA } },
	{ Rf,	"Ruther­fordium",	"Rf",	267,	true,  {	kNA,	kNA,	kNA,	157,	140,	131,	kNA } },
	{ Db,	"Dub­nium",			"Db",	268,	true,  {	kNA,	kNA,	kNA,	149,	136,	126,	kNA } },
	{ Sg,	"Sea­borgium",		"Sg",	269,	true,  {	kNA,	kNA,	kNA,	143,	128,	121,	kNA } },
	{ Bh,	"Bohr­ium",			"Bh",	270,	true,  {	kNA,	kNA,	kNA,	141,	128,	119,	kNA } },
	{ Hs,	"Has­sium",			"Hs",	277,	true,  {	kNA,	kNA,	kNA,	134,	125,	118,	kNA } },
	{ Mt,	"Meit­nerium",		"Mt",	278,	true,  {	kNA,	kNA,	kNA,	129,	125,	113,	kNA } },
	{ Ds,	"Darm­stadtium",		"Ds",	281,	true,  {	kNA,	kNA,	kNA,	128,	116,	112,	kNA } },
	{ Rg,	"Roent­genium",		"Rg",	282,	true,  {	kNA,	kNA,	kNA,	121,	116,	118,	kNA } },
	{ Cn,	"Coper­nicium",		"Cn",	285,	true,  {	kNA,	kNA,	kNA,	122,	137,	130,	kNA } },
	{ Nh,	"Nihon­ium",			"Nh",	286,	true,  {	kNA,	kNA,	kNA,	136,	kNA,	kNA,	kNA } },
	{ Fl,	"Flerov­ium",		"Fl",	289,	true,  {	kNA,	kNA,	kNA,	143,	kNA,	kNA,	kNA } },
	{ Mc,	"Moscov­ium",		"Mc",	290,	true,  {	kNA,	kNA,	kNA,	162,	kNA,	kNA,	kNA } },
	{ Lv,	"Liver­morium",		"Lv",	293,	true,  {	kNA,	kNA,	kNA,	175,	kNA,	kNA,	kNA } },
	{ Ts,	"Tenness­ine",		"Ts",	294,	true,  {	kNA,	kNA,	kNA,	165,	kNA,	kNA,	kNA } },
	{ Og,	"Oga­nesson",		"Og",	294,	true,  {	kNA,	kNA,	kNA,	157,	kNA,	kNA,	kNA } },
	{ Ce,	"Cerium",			"Ce",	140.12,	true,  {	kNA,	185,	kNA,	163,	137,	131,	kNA } },
	{ Pr,	"Praseo­dymium",		"Pr",	140.91,	true,  {	247,	185,	kNA,	176,	138,	128,	kNA } },
	{ Nd,	"Neo­dymium",		"Nd",	144.24,	true,  {	206,	185,	kNA,	174,	137,	kNA,	kNA } },
	{ Pm,	"Prome­thium",		"Pm",	145,	true,  {	205,	185,	kNA,	173,	135,	kNA,	kNA } },
	{ Sm,	"Sama­rium",			"Sm",	150.36,	true,  {	238,	185,	kNA,	172,	134,	kNA,	kNA } },
	{ Eu,	"Europ­ium",			"Eu",	151.96,	true,  {	231,	185,	kNA,	168,	134,	kNA,	kNA } },
	{ Gd,	"Gadolin­ium",		"Gd",	157.25,	true,  {	233,	180,	kNA,	169,	135,	132,	kNA } },
	{ Tb,	"Ter­bium",			"Tb",	158.93,	true,  {	225,	175,	kNA,	168,	135,	kNA,	kNA } },
	{ Dy,	"Dyspro­sium",		"Dy",	162.50,	true,  {	228,	175,	kNA,	167,	133,	kNA,	kNA } },
	{ Ho,	"Hol­mium",			"Ho",	164.93,	true,  {	226,	175,	kNA,	166,	133,	kNA,	kNA } },
	{ Er,	"Erbium",			"Er",	167.26,	true,  {	226,	175,	kNA,	165,	133,	kNA,	kNA } },
	{ Tm,	"Thulium",			"Tm",	168.93,	true,  {	222,	175,	kNA,	164,	131,	kNA,	kNA } },
	{ Yb,	"Ytter­bium",		"Yb",	173.05,	true,  {	222,	175,	kNA,	170,	129,	kNA,	kNA } },
	{ Lu,	"Lute­tium",			"Lu",	174.97,	true,  {	217,	175,	160,	162,	131,	131,	kNA } },
	{ Th,	"Thor­ium",			"Th",	232.04,	true,  {	kNA,	180,	kNA,	175,	143,	136,	kNA } },
	{ Pa,	"Protac­tinium",		"Pa",	231.04,	true,  {	kNA,	180,	kNA,	169,	138,	129,	kNA } },
	{ U,	"Ura­nium",			"U",	238.03,	true,  {	kNA,	175,	kNA,	170,	134,	118,	186 } },
	{ Np,	"Neptu­nium",		"Np",	237,	true,  {	kNA,	175,	kNA,	171,	136,	116,	kNA } },
	{ Pu,	"Pluto­nium",		"Pu",	244,	true,  {	kNA,	175,	kNA,	172,	135,	kNA,	kNA } },
	{ Am,	"Ameri­cium",		"Am",	243,	true,  {	kNA,	175,	kNA,	166,	135,	kNA,	kNA } },
	{ Cm,	"Curium",			"Cm",	247,	true,  {	kNA,	kNA,	kNA,	166,	136,	kNA,	kNA } },
	{ Bk,	"Berkel­ium",		"Bk",	247,	true,  {	kNA,	kNA,	kNA,	168,	139,	kNA,	kNA } },
	{ Cf,	"Califor­nium",		"Cf",	251,	true,  {	kNA,	kNA,	kNA,	168,	140,	kNA,	kNA } },
	{ Es,	"Einstei­nium",		"Es",	252,	true,  {	kNA,	kNA,	kNA,	165,	140,	kNA,	kNA } },
	{ Fm,	"Fer­mium",			"Fm",	257,	true,  {	kNA,	kNA,	kNA,	167,	kNA,	kNA,	kNA } },
	{ Md,	"Mende­levium",		"Md",	258,	true,  {	kNA,	kNA,	kNA,	173,	139,	kNA,	kNA } },
	{ No,	"Nobel­ium",			"No",	259,	true,  {	kNA,	kNA,	kNA,	176,	kNA,	kNA,	kNA } },
	{ Lr,	"Lawren­cium",		"Lr",	266,	true,  {	kNA,	kNA,	kNA,	161,	141,	kNA,	kNA } }
};

uint32 kKnownAtomsCount = sizeof(kKnownAtoms) / sizeof(AtomTypeInfo);

// --------------------------------------------------------------------
// AtomTypeTraits

AtomTypeTraits::AtomTypeTraits(const string& symbol)
	: mInfo(nullptr)
{
	for (auto& i: kKnownAtoms)
	{
		if (cif::iequals(i.symbol, symbol))
		{
			mInfo = &i;
			break;
		}
	}
	
	if (mInfo == nullptr)
		throw invalid_argument("Not a known element: " + symbol);
}

AtomTypeTraits::AtomTypeTraits(AtomType t)
{
	if (t < H or t > Lr)
		throw invalid_argument("atomType out of range");
	mInfo = &kKnownAtoms[t];
}

bool AtomTypeTraits::isElement(const string& symbol)
{
	bool result = false;
	
	for (auto& i: kKnownAtoms)
	{
		if (cif::iequals(i.symbol, symbol))
		{
			result = true;
			break;
		}
	}
	
	return result;
}

bool AtomTypeTraits::isMetal(const string& symbol)
{
	bool result = false;
	
	for (auto& i: kKnownAtoms)
	{
		if (cif::iequals(i.symbol, symbol))
		{
			result = i.metal;
			break;
		}
	}
	
	return result;
}
	
}
