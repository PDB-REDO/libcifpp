/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <cmath>

#include "cif++.hpp"

namespace cif
{

namespace data
{

const atom_type_info kKnownAtoms[] =
{
	{ Nn,	"Unknown",			"Nn",	0,			false, {	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },  //	0	Nn	 Unknown        
	{ H,	"Hydrogen",			"H",	1.008f,		false, {	53,		25,		37,		32,		kNA,	kNA,	120 } },  //	1	H	 Hydro­gen         
	{ He,	"Helium",			"He",	4.0026f,	false, {	31,		kNA,	32,		46,		kNA,	kNA,	140 } },  //	2	He	 He­lium           
	{ Li,	"Lithium",			"Li",	6.94f,		true,  {	167,	145,	134,	133,	124,	kNA,	182 } },  //	3	Li	 Lith­ium          
	{ Be,	"Beryllium",		"Be",	9.0122f,	true,  {	112,	105,	90,		102,	90,		85,		kNA } },  //	4	Be	 Beryl­lium        
	{ B,	"Boron",			"B",	10.81f,		true,  {	87,		85,		82,		85,		78,		73,		kNA } },  //	5	B	 Boron              
	{ C,	"Carbon",			"C",	12.011f,	false, {	67,		70,		77,		75,		67,		60,		170 } },  //	6	C	 Carbon             
	{ N,	"Nitrogen",			"N",	14.007f,	false, {	56,		65,		75,		71,		60,		54,		155 } },  //	7	N	 Nitro­gen         
	{ O,	"Oxygen",			"O",	15.999f,	false, {	48,		60,		73,		63,		57,		53,		152 } },  //	8	O	 Oxy­gen           
	{ F,	"Fluorine",			"F",	18.998f,	false, {	42,		50,		71,		64,		59,		53,		147 } },  //	9	F	 Fluor­ine         
	{ Ne,	"Neon",				"Ne",	20.180f,	false, {	38,		kNA,	69,		67,		96,		kNA,	154 } },  //	10	Ne	 Neon               
	{ Na,	"Sodium",			"Na",	22.990f,	true,  {	190,	180,	154,	155,	160,	kNA,	227 } },  //	11	Na	 So­dium           
	{ Mg,	"Magnesium",		"Mg",	24.305f,	true,  {	145,	150,	130,	139,	132,	127,	173 } },  //	12	Mg	 Magne­sium        
	{ Al,	"Aluminium",		"Al",	26.982f,	true,  {	118,	125,	118,	126,	113,	111,	kNA } },  //	13	Al	 Alumin­ium        
	{ Si,	"Silicon",			"Si",	28.085f,	true,  {	111,	110,	111,	116,	107,	102,	210 } },  //	14	Si	 Sili­con          
	{ P,	"Phosphorus",		"P",	30.974f,	false, {	98,		100,	106,	111,	102,	94,		180 } },  //	15	P	 Phos­phorus       
	{ S,	"Sulfur",			"S",	32.06f,		false, {	88,		100,	102,	103,	94,		95,		180 } },  //	16	S	 Sulfur             
	{ Cl,	"Chlorine",			"Cl",	35.45f,		false, {	79,		100,	99,		99,		95,		93,		175 } },  //	17	Cl	 Chlor­ine         
	{ Ar,	"Argon",			"Ar",	39.948f,	false, {	71,		kNA,	97,		96,		107,	96,		188 } },  //	18	Ar	 Argon              
	{ K,	"Potassium",		"K",	39.098f,	true,  {	243,	220,	196,	196,	193,	kNA,	275 } },  //	19	K	 Potas­sium        
	{ Ca,	"Calcium",			"Ca",	40.078f,	true,  {	194,	180,	174,	171,	147,	133,	kNA } },  //	20	Ca	 Cal­cium          
	{ Sc,	"Scandium",			"Sc",	44.956f,	true,  {	184,	160,	144,	148,	116,	114,	kNA } },  //	21	Sc	 Scan­dium         
	{ Ti,	"Titanium",			"Ti",	47.867f,	true,  {	176,	140,	136,	136,	117,	108,	kNA } },  //	22	Ti	 Tita­nium         
	{ V,	"Vanadium",			"V",	50.942f,	true,  {	171,	135,	125,	134,	112,	106,	kNA } },  //	23	V	 Vana­dium         
	{ Cr,	"Chromium",			"Cr",	51.996f,	true,  {	166,	140,	127,	122,	111,	103,	kNA } },  //	24	Cr	 Chrom­ium         
	{ Mn,	"Manganese",		"Mn",	54.938f,	true,  {	161,	140,	139,	119,	105,	103,	kNA } },  //	25	Mn	 Manga­nese        
	{ Fe,	"Iron",				"Fe",	55.845f,	true,  {	156,	140,	125,	116,	109,	102,	kNA } },  //	26	Fe	 Iron               
	{ Co,	"Cobalt",			"Co",	58.933f,	true,  {	152,	135,	126,	111,	103,	96,		kNA } },  //	27	Co	 Cobalt             
	{ Ni,	"Nickel",			"Ni",	58.693f,	true,  {	149,	135,	121,	110,	101,	101,	163 } },  //	28	Ni	 Nickel             
	{ Cu,	"Copper",			"Cu",	63.546f,	true,  {	145,	135,	138,	112,	115,	120,	140 } },  //	29	Cu	 Copper             
	{ Zn,	"Zinc",				"Zn",	65.38f,		true,  {	142,	135,	131,	118,	120,	kNA,	139 } },  //	30	Zn	 Zinc               
	{ Ga,	"Gallium",			"Ga",	69.723f,	true,  {	136,	130,	126,	124,	117,	121,	187 } },  //	31	Ga	 Gallium            
	{ Ge,	"Germanium",		"Ge",	72.630f,	true,  {	125,	125,	122,	121,	111,	114,	kNA } },  //	32	Ge	 Germa­nium        
	{ As,	"Arsenic",			"As",	74.922f,	true,  {	114,	115,	119,	121,	114,	106,	185 } },  //	33	As	 Arsenic            
	{ Se,	"Selenium",			"Se",	78.971f,	false, {	103,	115,	116,	116,	107,	107,	190 } },  //	34	Se	 Sele­nium         
	{ Br,	"Bromine",			"Br",	79.904f,	false, {	94,		115,	114,	114,	109,	110,	185 } },  //	35	Br	 Bromine            
	{ Kr,	"Krypton",			"Kr",	83.798f,	false, {	88,		kNA,	110,	117,	121,	108,	202 } },  //	36	Kr	 Kryp­ton          
	{ Rb,	"Rubidium",			"Rb",	85.468f,	true,  {	265,	235,	211,	210,	202,	kNA,	kNA } },  //	37	Rb	 Rubid­ium         
	{ Sr,	"Strontium",		"Sr",	87.62f,		true,  {	219,	200,	192,	185,	157,	139,	kNA } },  //	38	Sr	 Stront­ium        
	{ Y,	"Yttrium",			"Y",	88.906f,	true,  {	212,	180,	162,	163,	130,	124,	kNA } },  //	39	Y	 Yttrium            
	{ Zr,	"Zirconium",		"Zr",	91.224f,	true,  {	206,	155,	148,	154,	127,	121,	kNA } },  //	40	Zr	 Zirco­nium        
	{ Nb,	"Niobium",			"Nb",	92.906f,	true,  {	198,	145,	137,	147,	125,	116,	kNA } },  //	41	Nb	 Nio­bium          
	{ Mo,	"Molybdenum",		"Mo",	95.95f,		true,  {	190,	145,	145,	138,	121,	113,	kNA } },  //	42	Mo	 Molyb­denum       
	{ Tc,	"Technetium",		"Tc",	98,			true,  {	183,	135,	156,	128,	120,	110,	kNA } },  //	43	Tc	 Tech­netium       
	{ Ru,	"Ruthenium",		"Ru",	101.07f,	true,  {	178,	130,	126,	125,	114,	103,	kNA } },  //	44	Ru	 Ruthe­nium        
	{ Rh,	"Rhodium",			"Rh",	102.91f,	true,  {	173,	135,	135,	125,	110,	106,	kNA } },  //	45	Rh	 Rho­dium          
	{ Pd,	"Palladium",		"Pd",	106.42f,	true,  {	169,	140,	131,	120,	117,	112,	163 } },  //	46	Pd	 Pallad­ium        
	{ Ag,	"Silver",			"Ag",	107.87f,	true,  {	165,	160,	153,	128,	139,	137,	172 } },  //	47	Ag	 Silver             
	{ Cd,	"Cadmium",			"Cd",	112.41f,	true,  {	161,	155,	148,	136,	144,	kNA,	158 } },  //	48	Cd	 Cad­mium          
	{ In,	"Indium",			"In",	114.82f,	true,  {	156,	155,	144,	142,	136,	146,	193 } },  //	49	In	 Indium             
	{ Sn,	"Tin",				"Sn",	118.71f,	true,  {	145,	145,	141,	140,	130,	132,	217 } },  //	50	Sn	 Tin                
	{ Sb,	"Antimony",			"Sb",	121.76f,	false, {	133,	145,	138,	140,	133,	127,	kNA } },  //	51	Sb	 Anti­mony         
	{ Te,	"Tellurium",		"Te",	127.60f,	false, {	123,	140,	135,	136,	128,	121,	206 } },  //	52	Te	 Tellurium          
	{ I,	"Iodine",			"I",	126.90f,	false, {	115,	140,	133,	133,	129,	125,	198 } },  //	53	I	 Iodine             
	{ Xe,	"Xenon",			"Xe",	131.29f,	false, {	108,	kNA,	130,	131,	135,	122,	216 } },  //	54	Xe	 Xenon              
	{ Cs,	"Caesium",			"Cs",	132.91f,	true,  {	298,	260,	225,	232,	209,	kNA,	kNA } },  //	55	Cs	 Cae­sium          
	{ Ba,	"Barium",			"Ba",	137.33f,	true,  {	253,	215,	198,	196,	161,	149,	kNA } },  //	56	Ba	 Ba­rium           
	{ La,	"Lanthanum",		"La",	138.91f,	true,  {	kNA,	195,	169,	180,	139,	139,	kNA } },  //	57	La	 Lan­thanum        
	{ Ce,	"Cerium",			"Ce",	140.12f,	true,  {	kNA,	185,	kNA,	163,	137,	131,	kNA } },  //	58	Ce	 Cerium             
	{ Pr,	"Praseodymium",		"Pr",	140.91f,	true,  {	247,	185,	kNA,	176,	138,	128,	kNA } },  //	59	Pr	 Praseo­dymium     
	{ Nd,	"Neodymium",		"Nd",	144.24f,	true,  {	206,	185,	kNA,	174,	137,	kNA,	kNA } },  //	60	Nd	 Neo­dymium        
	{ Pm,	"Promethium",		"Pm",	145,		true,  {	205,	185,	kNA,	173,	135,	kNA,	kNA } },  //	61	Pm	 Prome­thium       
	{ Sm,	"Samarium",			"Sm",	150.36f,	true,  {	238,	185,	kNA,	172,	134,	kNA,	kNA } },  //	62	Sm	 Sama­rium         
	{ Eu,	"Europium",			"Eu",	151.96f,	true,  {	231,	185,	kNA,	168,	134,	kNA,	kNA } },  //	63	Eu	 Europ­ium         
	{ Gd,	"Gadolinium",		"Gd",	157.25f,	true,  {	233,	180,	kNA,	169,	135,	132,	kNA } },  //	64	Gd	 Gadolin­ium       
	{ Tb,	"Terbium",			"Tb",	158.93f,	true,  {	225,	175,	kNA,	168,	135,	kNA,	kNA } },  //	65	Tb	 Ter­bium          
	{ Dy,	"Dysprosium",		"Dy",	162.50f,	true,  {	228,	175,	kNA,	167,	133,	kNA,	kNA } },  //	66	Dy	 Dyspro­sium       
	{ Ho,	"Holmium",			"Ho",	164.93f,	true,  {	226,	175,	kNA,	166,	133,	kNA,	kNA } },  //	67	Ho	 Hol­mium          
	{ Er,	"Erbium",			"Er",	167.26f,	true,  {	226,	175,	kNA,	165,	133,	kNA,	kNA } },  //	68	Er	 Erbium             
	{ Tm,	"Thulium",			"Tm",	168.93f,	true,  {	222,	175,	kNA,	164,	131,	kNA,	kNA } },  //	69	Tm	 Thulium            
	{ Yb,	"Ytterbium",		"Yb",	173.05f,	true,  {	222,	175,	kNA,	170,	129,	kNA,	kNA } },  //	70	Yb	 Ytter­bium        
	{ Lu,	"Lutetium",			"Lu",	174.97f,	true,  {	217,	175,	160,	162,	131,	131,	kNA } },  //	71	Lu	 Lute­tium         
	{ Hf,	"Hafnium",			"Hf",	178.49f,	true,  {	208,	155,	150,	152,	128,	122,	kNA } },  //	72	Hf	 Haf­nium          
	{ Ta,	"Tantalum",			"Ta",	180.95f,	true,  {	200,	145,	138,	146,	126,	119,	kNA } },  //	73	Ta	 Tanta­lum         
	{ W,	"Tungsten",			"W",	183.84f,	true,  {	193,	135,	146,	137,	120,	115,	kNA } },  //	74	W	 Tung­sten         
	{ Re,	"Rhenium",			"Re",	186.21f,	true,  {	188,	135,	159,	131,	119,	110,	kNA } },  //	75	Re	 Rhe­nium          
	{ Os,	"Osmium",			"Os",	190.23f,	true,  {	185,	130,	128,	129,	116,	109,	kNA } },  //	76	Os	 Os­mium           
	{ Ir,	"Iridium",			"Ir",	192.22f,	true,  {	180,	135,	137,	122,	115,	107,	kNA } },  //	77	Ir	 Iridium            
	{ Pt,	"Platinum",			"Pt",	195.08f,	true,  {	177,	135,	128,	123,	112,	110,	175 } },  //	78	Pt	 Plat­inum         
	{ Au,	"Gold",				"Au",	196.97f,	true,  {	174,	135,	144,	124,	121,	123,	166 } },  //	79	Au	 Gold               
	{ Hg,	"Mercury",			"Hg",	200.59f,	true,  {	171,	150,	149,	133,	142,	kNA,	155 } },  //	80	Hg	 Mer­cury          
	{ Tl,	"Thallium",			"Tl",	204.38f,	true,  {	156,	190,	148,	144,	142,	150,	196 } },  //	81	Tl	 Thallium           
	{ Pb,	"Lead",				"Pb",	207.2f,		true,  {	154,	180,	147,	144,	135,	137,	202 } },  //	82	Pb	 Lead               
	{ Bi,	"Bismuth",			"Bi",	208.98f,	true,  {	143,	160,	146,	151,	141,	135,	kNA } },  //	83	Bi	 Bis­muth          
	{ Po,	"Polonium",			"Po",	209,		true,  {	135,	190,	kNA,	145,	135,	129,	kNA } },  //	84	Po	 Polo­nium         
	{ At,	"Astatine",			"At",	210,		false, {	127,	kNA,	kNA,	147,	138,	138,	kNA } },  //	85	At	 Asta­tine         
	{ Rn,	"Radon",			"Rn",	222,		false, {	120,	kNA,	145,	142,	145,	133,	kNA } },  //	86	Rn	 Radon              
	{ Fr,	"Francium",			"Fr",	223,		true,  {	kNA,	kNA,	kNA,	223,	218,	kNA,	kNA } },  //	87	Fr	 Fran­cium         
	{ Ra,	"Radium",			"Ra",	226,		true,  {	kNA,	215,	kNA,	201,	173,	159,	kNA } },  //	88	Ra	 Ra­dium           
	{ Ac,	"Actinium",			"Ac",	227,		true,  {	kNA,	195,	kNA,	186,	153,	140,	kNA } },  //	89	Ac	 Actin­ium         
	{ Th,	"Thorium",			"Th",	232.04f,	true,  {	kNA,	180,	kNA,	175,	143,	136,	kNA } },  //	90	Th	 Thor­ium          
	{ Pa,	"Protactinium",		"Pa",	231.04f,	true,  {	kNA,	180,	kNA,	169,	138,	129,	kNA } },  //	91	Pa	 Protac­tinium     
	{ U,	"Uranium",			"U",	238.03f,	true,  {	kNA,	175,	kNA,	170,	134,	118,	186 } },  //	92	U	 Ura­nium          
	{ Np,	"Neptunium",		"Np",	237,		true,  {	kNA,	175,	kNA,	171,	136,	116,	kNA } },  //	93	Np	 Neptu­nium        
	{ Pu,	"Plutonium",		"Pu",	244,		true,  {	kNA,	175,	kNA,	172,	135,	kNA,	kNA } },  //	94	Pu	 Pluto­nium        
	{ Am,	"Americium",		"Am",	243,		true,  {	kNA,	175,	kNA,	166,	135,	kNA,	kNA } },  //	95	Am	 Ameri­cium        
	{ Cm,	"Curium",			"Cm",	247,		true,  {	kNA,	kNA,	kNA,	166,	136,	kNA,	kNA } },  //	96	Cm	 Curium             
	{ Bk,	"Berkelium",		"Bk",	247,		true,  {	kNA,	kNA,	kNA,	168,	139,	kNA,	kNA } },  //	97	Bk	 Berkel­ium        
	{ Cf,	"Californium",		"Cf",	251,		true,  {	kNA,	kNA,	kNA,	168,	140,	kNA,	kNA } },  //	98	Cf	 Califor­nium      
	{ Es,	"Einsteinium",		"Es",	252,		true,  {	kNA,	kNA,	kNA,	165,	140,	kNA,	kNA } },  //	99	Es	 Einstei­nium      
	{ Fm,	"Fermium",			"Fm",	257,		true,  {	kNA,	kNA,	kNA,	167,	kNA,	kNA,	kNA } },  //	100	Fm	 Fer­mium          
	{ Md,	"Mendelevium",		"Md",	258,		true,  {	kNA,	kNA,	kNA,	173,	139,	kNA,	kNA } },  //	101	Md	 Mende­levium      
	{ No,	"Nobelium",			"No",	259,		true,  {	kNA,	kNA,	kNA,	176,	kNA,	kNA,	kNA } },  //	102	No	 Nobel­ium         
	{ Lr,	"Lawrencium",		"Lr",	266,		true,  {	kNA,	kNA,	kNA,	161,	141,	kNA,	kNA } },  //	103	Lr	 Lawren­cium       
	{ Rf,	"Rutherfordium",	"Rf",	267,		true,  {	kNA,	kNA,	kNA,	157,	140,	131,	kNA } },  //	104	Rf	 Ruther­fordium    
	{ Db,	"Dubnium",			"Db",	268,		true,  {	kNA,	kNA,	kNA,	149,	136,	126,	kNA } },  //	105	Db	 Dub­nium          
	{ Sg,	"Seaborgium",		"Sg",	269,		true,  {	kNA,	kNA,	kNA,	143,	128,	121,	kNA } },  //	106	Sg	 Sea­borgium       
	{ Bh,	"Bohrium",			"Bh",	270,		true,  {	kNA,	kNA,	kNA,	141,	128,	119,	kNA } },  //	107	Bh	 Bohr­ium          
	{ Hs,	"Hassium",			"Hs",	277,		true,  {	kNA,	kNA,	kNA,	134,	125,	118,	kNA } },  //	108	Hs	 Has­sium          
	{ Mt,	"Meitnerium",		"Mt",	278,		true,  {	kNA,	kNA,	kNA,	129,	125,	113,	kNA } },  //	109	Mt	 Meit­nerium       
	{ Ds,	"Darmstadtium",		"Ds",	281,		true,  {	kNA,	kNA,	kNA,	128,	116,	112,	kNA } },  //	110	Ds	 Darm­stadtium     
	{ Rg,	"Roentgenium",		"Rg",	282,		true,  {	kNA,	kNA,	kNA,	121,	116,	118,	kNA } },  //	111	Rg	 Roent­genium      
	{ Cn,	"Copernicium",		"Cn",	285,		true,  {	kNA,	kNA,	kNA,	122,	137,	130,	kNA } },  //	112	Cn	 Coper­nicium      
	{ Nh,	"Nihonium",			"Nh",	286,		true,  {	kNA,	kNA,	kNA,	136,	kNA,	kNA,	kNA } },  //	113	Nh	 Nihon­ium         
	{ Fl,	"Flerovium",		"Fl",	289,		true,  {	kNA,	kNA,	kNA,	143,	kNA,	kNA,	kNA } },  //	114	Fl	 Flerov­ium        
	{ Mc,	"Moscovium",		"Mc",	290,		true,  {	kNA,	kNA,	kNA,	162,	kNA,	kNA,	kNA } },  //	115	Mc	 Moscov­ium        
	{ Lv,	"Livermorium",		"Lv",	293,		true,  {	kNA,	kNA,	kNA,	175,	kNA,	kNA,	kNA } },  //	116	Lv	 Liver­morium      
	{ Ts,	"Tennessine",		"Ts",	294,		true,  {	kNA,	kNA,	kNA,	165,	kNA,	kNA,	kNA } },  //	117	Ts	 Tenness­ine       
	{ Og,	"Oganesson",		"Og",	294,		true,  {	kNA,	kNA,	kNA,	157,	kNA,	kNA,	kNA } },  //	118	Og	 Oga­nesson        

	{ D,	"Deuterium",		"D",	2.014f,		false, {	53,		25,		37,		32,		kNA,	kNA,	120 } },  //	1	D	 Deuterium         
};                                                                                                                                                    

uint32_t kKnownAtomsCount = sizeof(kKnownAtoms) / sizeof(atom_type_info);

// --------------------------------------------------------------------
// Crystal ionic radii, as taken from Wikipedia (https://en.m.wikipedia.org/wiki/Ionic_radius)

const struct ionic_radii
{
	atom_type type;
	float radii[11];
} kCrystalIonicRadii[] = {
	{ H,	{ kNA,	kNA,	208,	-4,		kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Hydrogen 
	{ Li,	{ kNA,	kNA,	kNA,	90,		kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Lithium 
	{ Be,	{ kNA,	kNA,	kNA,	kNA,	59,		kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Beryllium 
	{ B,	{ kNA,	kNA,	kNA,	kNA,	kNA,	41,		kNA,	kNA,	kNA,	kNA,	kNA } },	// Boron 
	{ C,	{ kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	30,		kNA,	kNA,	kNA,	kNA } },	// Carbon 
	{ N,	{ 132,	kNA,	kNA,	kNA,	kNA,	30,		kNA,	27,		kNA,	kNA,	kNA } },	// Nitrogen 
	{ O,	{ kNA,	126,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Oxygen 
	{ F,	{ kNA,	kNA,	119,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	22,		kNA } },	// Fluorine 
	{ Na,	{ kNA,	kNA,	kNA,	116,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Sodium 
	{ Mg,	{ kNA,	kNA,	kNA,	kNA,	86,		kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Magnesium 
	{ Al,	{ kNA,	kNA,	kNA,	kNA,	kNA,	67.5f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Aluminium 
	{ Si,	{ kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	54,		kNA,	kNA,	kNA,	kNA } },	// Silicon 
	{ P,	{ kNA,	kNA,	kNA,	kNA,	kNA,	58,		kNA,	52,		kNA,	kNA,	kNA } },	// Phosphorus 
	{ S,	{ kNA,	170,	kNA,	kNA,	kNA,	kNA,	51,		kNA,	43,		kNA,	kNA } },	// Sulfur 
	{ Cl,	{ kNA,	kNA,	181,	kNA,	kNA,	kNA,	kNA,	26,		kNA,	41,		kNA } },	// Chlorine 
	{ K,	{ kNA,	kNA,	kNA,	152,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Potassium 
	{ Ca,	{ kNA,	kNA,	kNA,	kNA,	114,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Calcium 
	{ Sc,	{ kNA,	kNA,	kNA,	kNA,	kNA,	88.5f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Scandium 
	{ Ti,	{ kNA,	kNA,	kNA,	kNA,	100,	81,		74.5f,	kNA,	kNA,	kNA,	kNA } },	// Titanium 
	{ V,	{ kNA,	kNA,	kNA,	kNA,	93,		78,		72,		68,		kNA,	kNA,	kNA } },	// Vanadium 
	{ Cr,	{ kNA,	kNA,	kNA,	kNA,	87,		75.5f,	69,		63,		58,		kNA,	kNA } },	// Chromium ls 
	// { Cr,{ kNA,	kNA,	kNA,	kNA,	94,		kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Chromium hs 
	{ Mn,	{ kNA,	kNA,	kNA,	kNA,	81,		72,		67,		47,		39.5f,	60,		kNA } },	// Manganese ls 
	// { Mn,{ kNA,	kNA,	kNA,	kNA,	97,		78.5f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Manganese hs 
	{ Fe,	{ kNA,	kNA,	kNA,	kNA,	75,		69,		72.5f,	kNA,	39,		kNA,	kNA } },	// Iron ls 
	// { Fe,{ kNA,	kNA,	kNA,	kNA,	92,		78.5f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Iron hs 
	{ Co,	{ kNA,	kNA,	kNA,	kNA,	79,		68.5f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Cobalt ls 
	// { Co,{ kNA,	kNA,	kNA,	kNA,	88.5f,	75,		67,		kNA,	kNA,	kNA,	kNA } },	// Cobalt hs 
	{ Ni,	{ kNA,	kNA,	kNA,	kNA,	83,		70,		62,		kNA,	kNA,	kNA,	kNA } },	// Nickel ls 
	// { Ni,{ kNA,	kNA,	kNA,	kNA,	kNA,	74,		kNA,	kNA,	kNA,	kNA,	kNA } },	// Nickel hs 
	{ Cu,	{ kNA,	kNA,	kNA,	91,		87,		68,		kNA,	kNA,	kNA,	kNA,	kNA } },	// Copper 
	{ Zn,	{ kNA,	kNA,	kNA,	kNA,	88	,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Zinc 
	{ Ga,	{ kNA,	kNA,	kNA,	kNA,	kNA,	76,		kNA,	kNA,	kNA,	kNA,	kNA } },	// Gallium 
	{ Ge,	{ kNA,	kNA,	kNA,	kNA,	87,		kNA,	67,		kNA,	kNA,	kNA,	kNA } },	// Germanium 
	{ As,	{ kNA,	kNA,	kNA,	kNA,	kNA,	72,		kNA,	60,		kNA,	kNA,	kNA } },	// Arsenic 
	{ Se,	{ kNA,	184,	kNA,	kNA,	kNA,	kNA,	64,		kNA,	56,		kNA,	kNA } },	// Selenium 
	{ Br,	{ kNA,	kNA,	182,	kNA,	kNA,	73,		kNA,	45,		kNA,	53,		kNA } },	// Bromine 
	{ Rb,	{ kNA,	kNA,	kNA,	166,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Rubidium 
	{ Sr,	{ kNA,	kNA,	kNA,	kNA,	132,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Strontium 
	{ Y,	{ kNA,	kNA,	kNA,	kNA,	kNA,	104,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Yttrium 
	{ Zr,	{ kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	86,		kNA,	kNA,	kNA,	kNA } },	// Zirconium 
	{ Nb,	{ kNA,	kNA,	kNA,	kNA,	kNA,	86,		82,		78,		kNA,	kNA,	kNA } },	// Niobium 
	{ Mo,	{ kNA,	kNA,	kNA,	kNA,	kNA,	83,		79,		75,		73,		kNA,	kNA } },	// Molybdenum 
	{ Tc,	{ kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	78.5f,	74,		kNA,	70,		kNA } },	// Technetium 
	{ Ru,	{ kNA,	kNA,	kNA,	kNA,	kNA,	82,		76,		70.5f,	kNA,	52,		150 } },	// Ruthenium 
	{ Rh,	{ kNA,	kNA,	kNA,	kNA,	kNA,	80.5f,	74,		69,		kNA,	kNA,	kNA } },	// Rhodium 
	{ Pd,	{ kNA,	kNA,	kNA,	73,		100,	90,		75.5f,	kNA,	kNA,	kNA,	kNA } },	// Palladium 
	{ Ag,	{ kNA,	kNA,	kNA,	129,	108,	89,		kNA,	kNA,	kNA,	kNA,	kNA } },	// Silver 
	{ Cd,	{ kNA,	kNA,	kNA,	kNA,	109,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Cadmium 
	{ In,	{ kNA,	kNA,	kNA,	kNA,	kNA,	94,		kNA,	kNA,	kNA,	kNA,	kNA } },	// Indium 
	{ Sn,	{ kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	83,		kNA,	kNA,	kNA,	kNA } },	// Tin 
	{ Sb,	{ kNA,	kNA,	kNA,	kNA,	kNA,	90,		kNA,	74,		kNA,	kNA,	kNA } },	// Antimony 
	{ Te,	{ kNA,	207,	kNA,	kNA,	kNA,	kNA,	111,	kNA,	70,		kNA,	kNA } },	// Tellurium 
	{ I,	{ kNA,	kNA,	206,	kNA,	kNA,	kNA,	kNA,	109,	kNA,	67,		kNA } },	// Iodine 
	{ Xe,	{ kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	62 } },		// Xenon 
	{ Cs,	{ kNA,	kNA,	kNA,	167,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Caesium 
	{ Ba,	{ kNA,	kNA,	kNA,	kNA,	149,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Barium 
	{ La,	{ kNA,	kNA,	kNA,	kNA,	kNA,	117.2f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Lanthanum 
	{ Ce,	{ kNA,	kNA,	kNA,	kNA,	kNA,	115,	101,	kNA,	kNA,	kNA,	kNA } },	// Cerium 
	{ Pr,	{ kNA,	kNA,	kNA,	kNA,	kNA,	113,	99,		kNA,	kNA,	kNA,	kNA } },	// Praseodymium 
	{ Nd,	{ kNA,	kNA,	kNA,	kNA,	143,	112.3f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Neodymium 
	{ Pm,	{ kNA,	kNA,	kNA,	kNA,	kNA,	111,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Promethium 
	{ Sm,	{ kNA,	kNA,	kNA,	kNA,	136,	109.8f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Samarium 
	{ Eu,	{ kNA,	kNA,	kNA,	kNA,	131,	108.7f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Europium 
	{ Gd,	{ kNA,	kNA,	kNA,	kNA,	kNA,	107.8f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Gadolinium 
	{ Tb,	{ kNA,	kNA,	kNA,	kNA,	kNA,	106.3f,	90,		kNA,	kNA,	kNA,	kNA } },	// Terbium 
	{ Dy,	{ kNA,	kNA,	kNA,	kNA,	121,	105.2f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Dysprosium 
	{ Ho,	{ kNA,	kNA,	kNA,	kNA,	kNA,	104.1f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Holmium 
	{ Er,	{ kNA,	kNA,	kNA,	kNA,	kNA,	103,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Erbium 
	{ Tm,	{ kNA,	kNA,	kNA,	kNA,	117,	102,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Thulium 
	{ Yb,	{ kNA,	kNA,	kNA,	kNA,	116,	100.8f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Ytterbium 
	{ Lu,	{ kNA,	kNA,	kNA,	kNA,	kNA,	100.1f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Lutetium 
	{ Hf,	{ kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	85,		kNA,	kNA,	kNA,	kNA } },	// Hafnium 
	{ Ta,	{ kNA,	kNA,	kNA,	kNA,	kNA,	86,		82,		78,		kNA,	kNA,	kNA } },	// Tantalum 
	{ W,	{ kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	80,		76,		74,		kNA,	kNA } },	// Tungsten 
	{ Re,	{ kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	77,		72,		69,		67,		kNA } },	// Rhenium 
	{ Os,	{ kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	77,		71.5f,	68.5f,	66.5f,	53 } },		// Osmium 
	{ Ir,	{ kNA,	kNA,	kNA,	kNA,	kNA,	82,		76.5f,	71,		kNA,	kNA,	kNA } },	// Iridium 
	{ Pt,	{ kNA,	kNA,	kNA,	kNA,	94,		kNA,	76.5f,	71,		kNA,	kNA,	kNA } },	// Platinum 
	{ Au,	{ kNA,	kNA,	kNA,	151,	kNA,	99,		kNA,	71,		kNA,	kNA,	kNA } },	// Gold 
	{ Hg,	{ kNA,	kNA,	kNA,	133,	116,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Mercury 
	{ Tl,	{ kNA,	kNA,	kNA,	164,	kNA,	102.5f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Thallium 
	{ Pb,	{ kNA,	kNA,	kNA,	kNA,	133,	kNA,	91.5f,	kNA,	kNA,	kNA,	kNA } },	// Lead 
	{ Bi,	{ kNA,	kNA,	kNA,	kNA,	kNA,	117,	kNA,	90,		kNA,	kNA,	kNA } },	// Bismuth 
	{ Po,	{ kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	108,	kNA,	81,		kNA,	kNA } },	// Polonium 
	{ At,	{ kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	76,		kNA } },	// Astatine 
	{ Fr,	{ kNA,	kNA,	kNA,	194,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Francium 
	{ Ra,	{ kNA,	kNA,	kNA,	kNA,	162,	kNA,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Radium 
	{ Ac,	{ kNA,	kNA,	kNA,	kNA,	kNA,	126,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Actinium 
	{ Th,	{ kNA,	kNA,	kNA,	kNA,	kNA,	kNA,	108,	kNA,	kNA,	kNA,	kNA } },	// Thorium 
	{ Pa,	{ kNA,	kNA,	kNA,	kNA,	kNA,	116,	104,	92,		kNA,	kNA,	kNA } },	// Protactinium 
	{ U,	{ kNA,	kNA,	kNA,	kNA,	kNA,	116.5f,	103,	90,		87,		kNA,	kNA } },	// Uranium 
	{ Np,	{ kNA,	kNA,	kNA,	kNA,	124,	115,	101,	89,		86,		85,		kNA } },	// Neptunium 
	{ Pu,	{ kNA,	kNA,	kNA,	kNA,	kNA,	114,	100,	88,		85,		kNA,	kNA } },	// Plutonium 
	{ Am,	{ kNA,	kNA,	kNA,	kNA,	140,	111.5f,	99,		kNA,	kNA,	kNA,	kNA } },	// Americium 
	{ Cm,	{ kNA,	kNA,	kNA,	kNA,	kNA,	111,	99,		kNA,	kNA,	kNA,	kNA } },	// Curium 
	{ Bk,	{ kNA,	kNA,	kNA,	kNA,	kNA,	110,	97,		kNA,	kNA,	kNA,	kNA } },	// Berkelium 
	{ Cf,	{ kNA,	kNA,	kNA,	kNA,	kNA,	109,	96.1f,	kNA,	kNA,	kNA,	kNA } },	// Californium 
	{ Es,	{ kNA,	kNA,	kNA,	kNA,	kNA,	92.8f,	kNA,	kNA,	kNA,	kNA,	kNA } },	// Einsteinium 	
}, kEffectiveIonicRadii[] = {
	{ H,	{ kNA, 	kNA, 	139.9f,	-18,	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Hydrogen 
	{ Li,	{ kNA, 	kNA, 	kNA, 	76,		kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Lithium 
	{ Be,	{ kNA, 	kNA, 	kNA, 	kNA, 	45,		kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Beryllium 
	{ B,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	27,		kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Boron 
	{ C,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	16,		kNA, 	kNA, 	kNA, 	kNA } }, 	// Carbon 
	{ N,	{ 146,	kNA, 	kNA, 	kNA, 	kNA, 	16,		kNA, 	13,		kNA, 	kNA, 	kNA } }, 	// Nitrogen 
	{ O,	{ kNA, 	140,	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Oxygen 
	{ F,	{ kNA, 	kNA, 	133,	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	8,		kNA } }, 	// Fluorine 
	{ Na,	{ kNA, 	kNA, 	kNA, 	102,	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Sodium 
	{ Mg,	{ kNA, 	kNA, 	kNA, 	kNA, 	72,		kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Magnesium 
	{ Al,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	53.5f,	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Aluminium 
	{ Si,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	40,		kNA, 	kNA, 	kNA, 	kNA } }, 	// Silicon 
	{ P,	{ 212,	kNA, 	kNA, 	kNA, 	kNA, 	44,		kNA, 	38,		kNA, 	kNA, 	kNA } }, 	// Phosphorus 
	{ S,	{ kNA, 	184,	kNA, 	kNA, 	kNA, 	kNA, 	37,		kNA, 	29,		kNA, 	kNA } }, 	// Sulfur 
	{ Cl,	{ kNA, 	kNA, 	181,	kNA, 	kNA, 	kNA, 	kNA, 	12,		kNA, 	27,		kNA } }, 	// Chlorine 
	{ K,	{ kNA, 	kNA, 	kNA, 	138,	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Potassium 
	{ Ca,	{ kNA, 	kNA, 	kNA, 	kNA, 	100,	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Calcium 
	{ Sc,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	74.5f,	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Scandium 
	{ Ti,	{ kNA, 	kNA, 	kNA, 	kNA, 	86,		67,		60.5f,	kNA, 	kNA, 	kNA, 	kNA } }, 	// Titanium 
	{ V,	{ kNA, 	kNA, 	kNA, 	kNA, 	79,		64,		58,		54,		kNA, 	kNA, 	kNA } }, 	// Vanadium 
	{ Cr,	{ kNA, 	kNA, 	kNA, 	kNA, 	73,		61.5f,	55,		49,		44,		kNA, 	kNA } }, 	// Chromium ls 
	{ Cr,	{ kNA, 	kNA, 	kNA, 	kNA, 	80,		kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Chromium hs 
	{ Mn,	{ kNA, 	kNA, 	kNA, 	kNA, 	67,		58,		53,		33,		25.5f,	46,		kNA } }, 	// Manganese ls 
	{ Mn,	{ kNA, 	kNA, 	kNA, 	kNA, 	83,		64.5f,	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Manganese hs 
	{ Fe,	{ kNA, 	kNA, 	kNA, 	kNA, 	61,		55,		58.5f,	kNA, 	25,		kNA, 	kNA } }, 	// Iron ls 
	{ Fe,	{ kNA, 	kNA, 	kNA, 	kNA, 	78,		64.5f,	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Iron hs 
	{ Co,	{ kNA, 	kNA, 	kNA, 	kNA, 	65,		54.5f,	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Cobalt ls 
	{ Co,	{ kNA, 	kNA, 	kNA, 	kNA, 	74.5f,	61,		53,		kNA, 	kNA, 	kNA, 	kNA } }, 	// Cobalt hs 
	{ Ni,	{ kNA, 	kNA, 	kNA, 	kNA, 	69,		56,		48,		kNA, 	kNA, 	kNA, 	kNA } }, 	// Nickel ls 
	{ Ni,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	60,		kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Nickel hs 
	{ Cu,	{ kNA, 	kNA, 	kNA, 	77,		73,		54, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Copper 
	{ Zn,	{ kNA, 	kNA, 	kNA, 	kNA, 	74,		kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Zinc 
	{ Ga,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	62,		kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Gallium 
	{ Ge,	{ kNA, 	kNA, 	kNA, 	kNA, 	73,		kNA, 	53,		kNA, 	kNA, 	kNA, 	kNA } }, 	// Germanium 
	{ As,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	58,		kNA, 	46,		kNA, 	kNA, 	kNA } }, 	// Arsenic 
	{ Se,	{ kNA, 	198,	kNA, 	kNA, 	kNA, 	kNA, 	50,		kNA, 	42,		kNA, 	kNA } }, 	// Selenium 
	{ Br,	{ kNA, 	kNA, 	196,	kNA, 	kNA, 	59,		kNA, 	31,		kNA, 	39,		kNA } }, 	// Bromine 
	{ Rb,	{ kNA, 	kNA, 	kNA, 	152,	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Rubidium 
	{ Sr,	{ kNA, 	kNA, 	kNA, 	kNA, 	118,	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Strontium 
	{ Y,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	90,		kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Yttrium 
	{ Zr,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	72,		kNA, 	kNA, 	kNA, 	kNA } }, 	// Zirconium 
	{ Nb,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	72,		68,		64,		kNA, 	kNA, 	kNA } }, 	// Niobium 
	{ Mo,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	69,		65,		61,		59,		kNA, 	kNA } }, 	// Molybdenum 
	{ Tc,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	64.5f,	60,		kNA, 	56,		kNA } }, 	// Technetium 
	{ Ru,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	68,		62,		56.5f,	kNA, 	38,		36 } },		// Ruthenium 
	{ Rh,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	66.5f,	60,		55,		kNA, 	kNA, 	kNA } }, 	// Rhodium 
	{ Pd,	{ kNA, 	kNA, 	kNA, 	59,		86,		76,		61.5f,	kNA, 	kNA, 	kNA, 	kNA } }, 	// Palladium 
	{ Ag,	{ kNA, 	kNA, 	kNA, 	115,	94,		75,		kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Silver 
	{ Cd,	{ kNA, 	kNA, 	kNA, 	kNA, 	95,		kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Cadmium 
	{ In,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	80,		kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Indium 
	{ Sn,	{ kNA, 	kNA, 	kNA, 	kNA, 	118,	kNA, 	69,		kNA, 	kNA, 	kNA, 	kNA } }, 	// Tin 
	{ Sb,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	76,		kNA, 	60,		kNA, 	kNA, 	kNA } }, 	// Antimony 
	{ Te,	{ kNA, 	221,	kNA, 	kNA, 	kNA, 	kNA, 	97,		kNA, 	56,		kNA, 	kNA } }, 	// Tellurium 
	{ I,	{ kNA, 	kNA, 	220,	kNA, 	kNA, 	kNA, 	kNA, 	95,		kNA, 	53,		kNA } }, 	// Iodine 
	{ Xe,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	48 } },		// Xenon 
	{ Cs,	{ kNA, 	kNA, 	kNA, 	167,	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Caesium 
	{ Ba,	{ kNA, 	kNA, 	kNA, 	kNA, 	135,	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Barium 
	{ La,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	103.2f,	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Lanthanum 
	{ Ce,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	101,	87,		kNA, 	kNA, 	kNA, 	kNA } }, 	// Cerium 
	{ Pr,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	99,		85,		kNA, 	kNA, 	kNA, 	kNA } }, 	// Praseodymium 
	{ Nd,	{ kNA, 	kNA, 	kNA, 	kNA, 	129,	98.3f,	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Neodymium 
	{ Pm,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	97,		kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Promethium 
	{ Sm,	{ kNA, 	kNA, 	kNA, 	kNA, 	122,	95.8f,	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Samarium 
	{ Eu,	{ kNA, 	kNA, 	kNA, 	kNA, 	117,	94.7f,	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Europium 
	{ Gd,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	93.5f,	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Gadolinium 
	{ Tb,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	92.3f,	76,		kNA, 	kNA, 	kNA, 	kNA } }, 	// Terbium 
	{ Dy,	{ kNA, 	kNA, 	kNA, 	kNA, 	107,	91.2f,	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Dysprosium 
	{ Ho,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	90.1f,	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Holmium 
	{ Er,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	89,		kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Erbium 
	{ Tm,	{ kNA, 	kNA, 	kNA, 	kNA, 	103,	88,		kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Thulium 
	{ Yb,	{ kNA, 	kNA, 	kNA, 	kNA, 	102,	86.8f,	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Ytterbium 
	{ Lu,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	86.1f,	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Lutetium 
	{ Hf,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	71,		kNA, 	kNA, 	kNA, 	kNA } }, 	// Hafnium 
	{ Ta,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	72,		68,		64,		kNA, 	kNA, 	kNA } }, 	// Tantalum 
	{ W,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	66,		62,		60,		kNA, 	kNA } }, 	// Tungsten 
	{ Re,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	63,		58,		55,		53,		kNA } }, 	// Rhenium 
	{ Os,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	63,		57.5f,	54.5f,	52.5f,	39 } },		// Osmium 
	{ Ir,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	68,		62.5f,	57,		kNA, 	kNA, 	kNA } }, 	// Iridium 
	{ Pt,	{ kNA, 	kNA, 	kNA, 	kNA, 	80,		kNA, 	62.5f,	57,		kNA, 	kNA, 	kNA } }, 	// Platinum 
	{ Au,	{ kNA, 	kNA, 	kNA, 	137,	kNA, 	85,		kNA, 	57,		kNA, 	kNA, 	kNA } }, 	// Gold 
	{ Hg,	{ kNA, 	kNA, 	kNA, 	119,	102,	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Mercury 
	{ Tl,	{ kNA, 	kNA, 	kNA, 	150,	kNA, 	88.5f,	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Thallium 
	{ Pb,	{ kNA, 	kNA, 	kNA, 	kNA, 	119,	kNA, 	77.5f,	kNA, 	kNA, 	kNA, 	kNA } }, 	// Lead 
	{ Bi,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	103,	kNA, 	76,		kNA, 	kNA, 	kNA } }, 	// Bismuth 
	{ Po,	{ kNA, 	223,	kNA, 	kNA, 	kNA, 	kNA, 	94,		kNA, 	67,		kNA, 	kNA } }, 	// Polonium 
	{ At,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	62,		kNA } }, 	// Astatine 
	{ Fr,	{ kNA, 	kNA, 	kNA, 	180,	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Francium 
	{ Ra,	{ kNA, 	kNA, 	kNA, 	kNA, 	148,	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } }, 	// Radium 
	{ Ac,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	106.5f, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA } },		// Actinium 
	{ Th,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	94,		kNA, 	kNA, 	kNA, 	kNA } }, 	// Thorium 
	{ Pa,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	104,	90,		78,		kNA, 	kNA, 	kNA } }, 	// Protactinium 
	{ U,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	102.5f,	89,		76,		73,		kNA, 	kNA } }, 	// Uranium 
	{ Np,	{ kNA, 	kNA, 	kNA, 	kNA, 	110,	101,	87,		75,		72,		71,		kNA } }, 	// Neptunium 
	{ Pu,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	100,	86,		74,		71,		kNA, 	kNA } }, 	// Plutonium 
	{ Am,	{ kNA, 	kNA, 	kNA, 	kNA, 	126,	97.5f,	85,		kNA, 	kNA, 	kNA, 	kNA } }, 	// Americium 
	{ Cm,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	97,		85,		kNA, 	kNA, 	kNA, 	kNA } }, 	// Curium 
	{ Bk,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	96,		83,		kNA, 	kNA, 	kNA, 	kNA } }, 	// Berkelium 
	{ Cf,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	95,		82.1f,	kNA, 	kNA, 	kNA, 	kNA } }, 	// Californium 
	{ Es,	{ kNA, 	kNA, 	kNA, 	kNA, 	kNA, 	83.5f,	kNA, 	kNA,	kNA,	kNA, 	kNA } }, 	// Einsteinium 	
};


// --------------------------------------------------------------------
// The coefficients from Waasmaier & Kirfel (1995), Acta Cryst. A51, 416-431.

struct SFDataArrayElement
{
	atom_type				symbol;
	int8_t					charge;
	atom_type_traits::SFData	sf;
};

SFDataArrayElement kWKSFData[] = {
	{ H,	0,	{{  0.413048,  0.294953,  0.187491,  0.080701,  0.023736,  0.000049},
				{ 15.569946, 32.398468,  5.711404, 61.889874,  1.334118,  0.000000}}},
	{ He,	0,	{{  0.732354,  0.753896,  0.283819,  0.190003,  0.039139,  0.000487},
				{ 11.553918,  4.595831,  1.546299, 26.463964,  0.377523,  0.000000}}},
	{ Li,	0,	{{  0.974637,  0.158472,  0.811855,  0.262416,  0.790108,  0.002542},
				{  4.334946,  0.342451, 97.102966,201.363831,  1.409234,  0.000000}}},
	{ Be,	0,	{{  1.533712,  0.638283,  0.601052,  0.106139,  1.118414,  0.002511},
				{ 42.662079,  0.595420, 99.106499,  0.151340,  1.843093,  0.000000}}},
	{ B,	0,	{{  2.085185,  1.064580,  1.062788,  0.140515,  0.641784,  0.003823},
				{ 23.494068,  1.137894, 61.238976,  0.114886,  0.399036,  0.000000}}},
	{ C,	0,	{{  2.657506,  1.078079,  1.490909, -4.241070,  0.713791,  4.297983},
				{ 14.780758,  0.776775, 42.086842, -0.000294,  0.239535,  0.000000}}},
	{ N,	0,	{{ 11.893780,  3.277479,  1.858092,  0.858927,  0.912985,-11.804902},
				{  0.000158, 10.232723, 30.344690,  0.656065,  0.217287,  0.000000}}},
	{ O,	0,	{{  2.960427,  2.508818,  0.637853,  0.722838,  1.142756,  0.027014},
				{ 14.182259,  5.936858,  0.112726, 34.958481,  0.390240,  0.000000}}},
	{ F,	0,	{{  3.511943,  2.772244,  0.678385,  0.915159,  1.089261,  0.032557},
				{ 10.687859,  4.380466,  0.093982, 27.255203,  0.313066,  0.000000}}},
	{ Ne,	0,	{{  4.183749,  2.905726,  0.520513,  1.135641,  1.228065,  0.025576},
				{  8.175457,  3.252536,  0.063295, 21.813910,  0.224952,  0.000000}}},
	{ Na,	0,	{{  4.910127,  3.081783,  1.262067,  1.098938,  0.560991,  0.079712},
				{  3.281434,  9.119178,  0.102763,132.013947,  0.405878,  0.000000}}},
	{ Mg,	0,	{{  4.708971,  1.194814,  1.558157,  1.170413,  3.239403,  0.126842},
				{  4.875207,108.506081,  0.111516, 48.292408,  1.928171,  0.000000}}},
	{ Al,	0,	{{  4.730796,  2.313951,  1.541980,  1.117564,  3.154754,  0.139509},
				{  3.628931, 43.051167,  0.095960,108.932388,  1.555918,  0.000000}}},
	{ Si,	0,	{{  5.275329,  3.191038,  1.511514,  1.356849,  2.519114,  0.145073},
				{  2.631338, 33.730728,  0.081119, 86.288643,  1.170087,  0.000000}}},
	{ P,	0,	{{  1.950541,  4.146930,  1.494560,  1.522042,  5.729711,  0.155233},
				{  0.908139, 27.044952,  0.071280, 67.520187,  1.981173,  0.000000}}},
	{ S,	0,	{{  6.372157,  5.154568,  1.473732,  1.635073,  1.209372,  0.154722},
				{  1.514347, 22.092527,  0.061373, 55.445175,  0.646925,  0.000000}}},
	{ Cl,	0,	{{  1.446071,  6.870609,  6.151801,  1.750347,  0.634168,  0.146773},
				{  0.052357,  1.193165, 18.343416, 46.398396,  0.401005,  0.000000}}},
	{ Ar,	0,	{{  7.188004,  6.638454,  0.454180,  1.929593,  1.523654,  0.265954},
				{  0.956221, 15.339877, 15.339862, 39.043823,  0.062409,  0.000000}}},
	{ K,	0,	{{  8.163991,  7.146945,  1.070140,  0.877316,  1.486434,  0.253614},
				{ 12.816323,  0.808945,210.327011, 39.597652,  0.052821,  0.000000}}},
	{ Ca,	0,	{{  8.593655,  1.477324,  1.436254,  1.182839,  7.113258,  0.196255},
				{ 10.460644,  0.041891, 81.390381,169.847839,  0.688098,  0.000000}}},
	{ Sc,	0,	{{  1.476566,  1.487278,  1.600187,  9.177463,  7.099750,  0.157765},
				{ 53.131023,  0.035325,137.319489,  9.098031,  0.602102,  0.000000}}},
	{ Ti,	0,	{{  9.818524,  1.522646,  1.703101,  1.768774,  7.082555,  0.102473},
				{  8.001879,  0.029763, 39.885422,120.157997,  0.532405,  0.000000}}},
	{ V,	0,	{{ 10.473575,  1.547881,  1.986381,  1.865616,  7.056250,  0.067744},
				{  7.081940,  0.026040, 31.909672,108.022842,  0.474882,  0.000000}}},
	{ Cr,	0,	{{ 11.007069,  1.555477,  2.985293,  1.347855,  7.034779,  0.065510},
				{  6.366281,  0.023987, 23.244839,105.774498,  0.429369,  0.000000}}},
	{ Mn,	0,	{{ 11.709542,  1.733414,  2.673141,  2.023368,  7.003180, -0.147293},
				{  5.597120,  0.017800, 21.788420, 89.517914,  0.383054,  0.000000}}},
	{ Fe,	0,	{{ 12.311098,  1.876623,  3.066177,  2.070451,  6.975185, -0.304931},
				{  5.009415,  0.014461, 18.743040, 82.767876,  0.346506,  0.000000}}},
	{ Co,	0,	{{ 12.914510,  2.481908,  3.466894,  2.106351,  6.960892, -0.936572},
				{  4.507138,  0.009126, 16.438129, 76.987320,  0.314418,  0.000000}}},
	{ Ni,	0,	{{ 13.521865,  6.947285,  3.866028,  2.135900,  4.284731, -2.762697},
				{  4.077277,  0.286763, 14.622634, 71.966080,  0.004437,  0.000000}}},
	{ Cu,	0,	{{ 14.014192,  4.784577,  5.056806,  1.457971,  6.932996, -3.254477},
				{  3.738280,  0.003744, 13.034982, 72.554794,  0.265666,  0.000000}}},
	{ Zn,	0,	{{ 14.741002,  6.907748,  4.642337,  2.191766, 38.424042,-36.915829},
				{  3.388232,  0.243315, 11.903689, 63.312130,  0.000397,  0.000000}}},
	{ Ga,	0,	{{ 15.758946,  6.841123,  4.121016,  2.714681,  2.395246, -0.847395},
				{  3.121754,  0.226057, 12.482196, 66.203621,  0.007238,  0.000000}}},
	{ Ge,	0,	{{ 16.540613,  1.567900,  3.727829,  3.345098,  6.785079,  0.018726},
				{  2.866618,  0.012198, 13.432163, 58.866047,  0.210974,  0.000000}}},
	{ As,	0,	{{ 17.025642,  4.503441,  3.715904,  3.937200,  6.790175, -2.984117},
				{  2.597739,  0.003012, 14.272119, 50.437996,  0.193015,  0.000000}}},
	{ Se,	0,	{{ 17.354071,  4.653248,  4.259489,  4.136455,  6.749163, -3.160982},
				{  2.349787,  0.002550, 15.579460, 45.181202,  0.177432,  0.000000}}},
	{ Br,	0,	{{ 17.550570,  5.411882,  3.937180,  3.880645,  6.707793, -2.492088},
				{  2.119226, 16.557184,  0.002481, 42.164009,  0.162121,  0.000000}}},
	{ Kr,	0,	{{ 17.655279,  6.848105,  4.171004,  3.446760,  6.685200, -2.810592},
				{  1.908231, 16.606236,  0.001598, 39.917473,  0.146896,  0.000000}}},
	{ Rb,	0,	{{  8.123134,  2.138042,  6.761702,  1.156051, 17.679546,  1.139548},
				{ 15.142385, 33.542667,  0.129372,224.132507,  1.713368,  0.000000}}},
	{ Sr,	0,	{{ 17.730219,  9.795867,  6.099763,  2.620025,  0.600053,  1.140251},
				{  1.563060, 14.310868,  0.120574,135.771317,  0.120574,  0.000000}}},
	{ Y,	0,	{{ 17.792040, 10.253252,  5.714949,  3.170516,  0.918251,  1.131787},
				{  1.429691, 13.132816,  0.112173,108.197029,  0.112173,  0.000000}}},
	{ Zr,	0,	{{ 17.859772, 10.911038,  5.821115,  3.512513,  0.746965,  1.124859},
				{  1.310692, 12.319285,  0.104353, 91.777542,  0.104353,  0.000000}}},
	{ Nb,	0,	{{ 17.958399, 12.063054,  5.007015,  3.287667,  1.531019,  1.123452},
				{  1.211590, 12.246687,  0.098615, 75.011948,  0.098615,  0.000000}}},
	{ Mo,	0,	{{  6.236218, 17.987711, 12.973127,  3.451426,  0.210899,  1.108770},
				{  0.090780,  1.108310, 11.468720, 66.684151,  0.090780,  0.000000}}},
	{ Tc,	0,	{{ 17.840963,  3.428236,  1.373012, 12.947364,  6.335469,  1.074784},
				{  1.005729, 41.901382,119.320541,  9.781542,  0.083391,  0.000000}}},
	{ Ru,	0,	{{  6.271624, 17.906738, 14.123269,  3.746008,  0.908235,  1.043992},
				{  0.077040,  0.928222,  9.555345, 35.860680,123.552246,  0.000000}}},
	{ Rh,	0,	{{  6.216648, 17.919739,  3.854252,  0.840326, 15.173498,  0.995452},
				{  0.070789,  0.856121, 33.889484,121.686691,  9.029517,  0.000000}}},
	{ Pd,	0,	{{  6.121511,  4.784063, 16.631683,  4.318258, 13.246773,  0.883099},
				{  0.062549,  0.784031,  8.751391, 34.489983,  0.784031,  0.000000}}},
	{ Ag,	0,	{{  6.073874, 17.155437,  4.173344,  0.852238, 17.988686,  0.756603},
				{  0.055333,  7.896512, 28.443739,110.376106,  0.716809,  0.000000}}},
	{ Cd,	0,	{{  6.080986, 18.019468,  4.018197,  1.303510, 17.974669,  0.603504},
				{  0.048990,  7.273646, 29.119284, 95.831207,  0.661231,  0.000000}}},
	{ In,	0,	{{  6.196477, 18.816183,  4.050479,  1.638929, 17.962912,  0.333097},
				{  0.042072,  6.695665, 31.009790,103.284348,  0.610714,  0.000000}}},
	{ Sn,	0,	{{ 19.325171,  6.281571,  4.498866,  1.856934, 17.917318,  0.119024},
				{  6.118104,  0.036915, 32.529045, 95.037186,  0.565651,  0.000000}}},
	{ Sb,	0,	{{  5.394956,  6.549570, 19.650681,  1.827820, 17.867832, -0.290506},
				{ 33.326523,  0.030974,  5.564929, 87.130966,  0.523992,  0.000000}}},
	{ Te,	0,	{{  6.660302,  6.940756, 19.847015,  1.557175, 17.802427, -0.806668},
				{ 33.031654,  0.025750,  5.065547, 84.101616,  0.487660,  0.000000}}},
	{ I,	0,	{{ 19.884502,  6.736593,  8.110516,  1.170953, 17.548716, -0.448811},
				{  4.628591,  0.027754, 31.849096, 84.406387,  0.463550,  0.000000}}},
	{ Xe,	0,	{{ 19.978920, 11.774945,  9.332182,  1.244749, 17.737501, -6.065902},
				{  4.143356,  0.010142, 28.796200, 75.280685,  0.413616,  0.000000}}},
	{ Cs,	0,	{{ 17.418674,  8.314444, 10.323193,  1.383834, 19.876251, -2.322802},
				{  0.399828,  0.016872, 25.605827,233.339676,  3.826915,  0.000000}}},
	{ Ba,	0,	{{ 19.747343, 17.368477, 10.465718,  2.592602, 11.003653, -5.183497},
				{  3.481823,  0.371224, 21.226641,173.834274,  0.010719,  0.000000}}},
	{ La,	0,	{{ 19.966019, 27.329655, 11.018425,  3.086696, 17.335455,-21.745489},
				{  3.197408,  0.003446, 19.955492,141.381973,  0.341817,  0.000000}}},
	{ Ce,	0,	{{ 17.355122, 43.988499, 20.546650,  3.130670, 11.353665,-38.386017},
				{  0.328369,  0.002047,  3.088196,134.907654, 18.832960,  0.000000}}},
	{ Pr,	0,	{{ 21.551311, 17.161730, 11.903859,  2.679103,  9.564197, -3.871068},
				{  2.995675,  0.312491, 17.716705,152.192825,  0.010468,  0.000000}}},
	{ Nd,	0,	{{ 17.331244, 62.783924, 12.160097,  2.663483, 22.239950,-57.189842},
				{  0.300269,  0.001320, 17.026001,148.748993,  2.910268,  0.000000}}},
	{ Pm,	0,	{{ 17.286388, 51.560162, 12.478557,  2.675515, 22.960947,-45.973682},
				{  0.286620,  0.001550, 16.223755,143.984512,  2.796480,  0.000000}}},
	{ Sm,	0,	{{ 23.700363, 23.072214, 12.777782,  2.684217, 17.204367,-17.452166},
				{  2.689539,  0.003491, 15.495437,139.862473,  0.274536,  0.000000}}},
	{ Eu,	0,	{{ 17.186195, 37.156837, 13.103387,  2.707246, 24.419271,-31.586687},
				{  0.261678,  0.001995, 14.787360,134.816299,  2.581883,  0.000000}}},
	{ Gd,	0,	{{ 24.898117, 17.104952, 13.222581,  3.266152, 48.995213,-43.505684},
				{  2.435028,  0.246961, 13.996325,110.863091,  0.001383,  0.000000}}},
	{ Tb,	0,	{{ 25.910013, 32.344139, 13.765117,  2.751404, 17.064405,-26.851971},
				{  2.373912,  0.002034, 13.481969,125.836510,  0.236916,  0.000000}}},
	{ Dy,	0,	{{ 26.671785, 88.687576, 14.065445,  2.768497, 17.067781,-83.279831},
				{  2.282593,  0.000665, 12.920230,121.937187,  0.225531,  0.000000}}},
	{ Ho,	0,	{{ 27.150190, 16.999819, 14.059334,  3.386979, 46.546471,-41.165253},
				{  2.169660,  0.215414, 12.213148,100.506783,  0.001211,  0.000000}}},
	{ Er,	0,	{{ 28.174887, 82.493271, 14.624002,  2.802756, 17.018515,-77.135223},
				{  2.120995,  0.000640, 11.915256,114.529938,  0.207519,  0.000000}}},
	{ Tm,	0,	{{ 28.925894, 76.173798, 14.904704,  2.814812, 16.998117,-70.839813},
				{  2.046203,  0.000656, 11.465375,111.411980,  0.199376,  0.000000}}},
	{ Yb,	0,	{{ 29.676760, 65.624069, 15.160854,  2.830288, 16.997850,-60.313812},
				{  1.977630,  0.000720, 11.044622,108.139153,  0.192110,  0.000000}}},
	{ Lu,	0,	{{ 30.122866, 15.099346, 56.314899,  3.540980, 16.943729,-51.049416},
				{  1.883090, 10.342764,  0.000780, 89.559250,  0.183849,  0.000000}}},
	{ Hf,	0,	{{ 30.617033, 15.145351, 54.933548,  4.096253, 16.896156,-49.719837},
				{  1.795613,  9.934469,  0.000739, 76.189705,  0.175914,  0.000000}}},
	{ Ta,	0,	{{ 31.066359, 15.341823, 49.278297,  4.577665, 16.828321,-44.119026},
				{  1.708732,  9.618455,  0.000760, 66.346199,  0.168002,  0.000000}}},
	{ W,	0,	{{ 31.507900, 15.682498, 37.960129,  4.885509, 16.792112,-32.864574},
				{  1.629485,  9.446448,  0.000898, 59.980675,  0.160798,  0.000000}}},
	{ Re,	0,	{{ 31.888456, 16.117104, 42.390297,  5.211669, 16.767591,-37.412682},
				{  1.549238,  9.233474,  0.000689, 54.516373,  0.152815,  0.000000}}},
	{ Os,	0,	{{ 32.210297, 16.678440, 48.559906,  5.455839, 16.735533,-43.677956},
				{  1.473531,  9.049695,  0.000519, 50.210201,  0.145771,  0.000000}}},
	{ Ir,	0,	{{ 32.004436,  1.975454, 17.070105, 15.939454,  5.990003,  4.018893},
				{  1.353767, 81.014175,  0.128093,  7.661196, 26.659403,  0.000000}}},
	{ Pt,	0,	{{ 31.273891, 18.445440, 17.063745,  5.555933,  1.575270,  4.050394},
				{  1.316992,  8.797154,  0.124741, 40.177994,  1.316997,  0.000000}}},
	{ Au,	0,	{{ 16.777390, 19.317156, 32.979683,  5.595453, 10.576854, -6.279078},
				{  0.122737,  8.621570,  1.256902, 38.008820,  0.000601,  0.000000}}},
	{ Hg,	0,	{{ 16.839890, 20.023823, 28.428564,  5.881564,  4.714706,  4.076478},
				{  0.115905,  8.256927,  1.195250, 39.247227,  1.195250,  0.000000}}},
	{ Tl,	0,	{{ 16.630795, 19.386616, 32.808571,  1.747191,  6.356862,  4.066939},
				{  0.110704,  7.181401,  1.119730, 90.660263, 26.014978,  0.000000}}},
	{ Pb,	0,	{{ 16.419567, 32.738590,  6.530247,  2.342742, 19.916475,  4.049824},
				{  0.105499,  1.055049, 25.025890, 80.906593,  6.664449,  0.000000}}},
	{ Bi,	0,	{{ 16.282274, 32.725136,  6.678302,  2.694750, 20.576559,  4.040914},
				{  0.101180,  1.002287, 25.714146, 77.057549,  6.291882,  0.000000}}},
	{ Po,	0,	{{ 16.289164, 32.807171, 21.095163,  2.505901,  7.254589,  4.046556},
				{  0.098121,  0.966265,  6.046622, 76.598068, 28.096128,  0.000000}}},
	{ At,	0,	{{ 16.011461, 32.615547,  8.113899,  2.884082, 21.377867,  3.995684},
				{  0.092639,  0.904416, 26.543257, 68.372963,  5.499512,  0.000000}}},
	{ Rn,	0,	{{ 16.070229, 32.641106, 21.489658,  2.299218,  9.480184,  4.020977},
				{  0.090437,  0.876409,  5.239687, 69.188477, 27.632641,  0.000000}}},
	{ Fr,	0,	{{ 16.007385, 32.663830, 21.594351,  1.598497, 11.121192,  4.003472},
				{  0.087031,  0.840187,  4.954467,199.805801, 26.905106,  0.000000}}},
	{ Ra,	0,	{{ 32.563690, 21.396671, 11.298093,  2.834688, 15.914965,  3.981773},
				{  0.801980,  4.590666, 22.758972,160.404388,  0.083544,  0.000000}}},
	{ Ac,	0,	{{ 15.914053, 32.535042, 21.553976, 11.433394,  3.612409,  3.939212},
				{  0.080511,  0.770669,  4.352206, 21.381622,130.500748,  0.000000}}},
	{ Th,	0,	{{ 15.784024, 32.454899, 21.849222,  4.239077, 11.736191,  3.922533},
				{  0.077067,  0.735137,  4.097976,109.464111, 20.512138,  0.000000}}},
	{ Pa,	0,	{{ 32.740208, 21.973675, 12.957398,  3.683832, 15.744058,  3.886066},
				{  0.709545,  4.050881, 19.231543,117.255005,  0.074040,  0.000000}}},
	{ U,	0,	{{ 15.679275, 32.824306, 13.660459,  3.687261, 22.279434,  3.854444},
				{  0.071206,  0.681177, 18.236156,112.500038,  3.930325,  0.000000}}},
	{ Np,	0,	{{ 32.999901, 22.638077, 14.219973,  3.672950, 15.683245,  3.769391},
				{  0.657086,  3.854918, 17.435474,109.464485,  0.068033,  0.000000}}},
	{ Pu,	0,	{{ 33.281178, 23.148544, 15.153755,  3.031492, 15.704215,  3.664200},
				{  0.634999,  3.856168, 16.849735,121.292038,  0.064857,  0.000000}}},
	{ Am,	0,	{{ 33.435162, 23.657259, 15.576339,  3.027023, 15.746100,  3.541160},
				{  0.612785,  3.792942, 16.195778,117.757004,  0.061755,  0.000000}}},
	{ Cm,	0,	{{ 15.804837, 33.480801, 24.150198,  3.655563, 15.499866,  3.390840},
				{  0.058619,  0.590160,  3.674720,100.736191, 15.408296,  0.000000}}},
	{ Bk,	0,	{{ 15.889072, 33.625286, 24.710381,  3.707139, 15.839268,  3.213169},
				{  0.055503,  0.569571,  3.615472, 97.694786, 14.754303,  0.000000}}},
	{ Cf,	0,	{{ 33.794075, 25.467693, 16.048487,  3.657525, 16.008982,  3.005326},
				{  0.550447,  3.581973, 14.357388, 96.064972,  0.052450,  0.000000}}},
	{ H,	-1,	{{  0.702260,  0.763666,  0.248678,  0.261323,  0.023017,  0.000425},
				{ 23.945604, 74.897919,  6.773289,233.583450,  1.337531,  0.000000}}},
	{ Li,	+1,	{{  0.432724,  0.549257,  0.376575, -0.336481,  0.976060,  0.001764},
				{  0.260367,  1.042836,  7.885294,  0.260368,  3.042539,  0.000000}}},
	{ Be,	+2,	{{  3.055430, -2.372617,  1.044914,  0.544233,  0.381737, -0.653773},
				{  0.001226,  0.001227,  1.542106,  0.456279,  4.047479,  0.000000}}},
	{ C, atom_type_traits::kWKSFVal,
				{{  1.258489,  0.728215,  1.119856,  2.168133,  0.705239,  0.019722},
				{ 10.683769,  0.208177,  0.836097, 24.603704, 58.954273,  0.000000}}},
	{ O,	-1,	{{  3.106934,  3.235142,  1.148886,  0.783981,  0.676953,  0.046136},
				{ 19.868080,  6.960252,  0.170043, 65.693512,  0.630757,  0.000000}}},
	{ O,	-2,	{{  3.990247,  2.300563,  0.607200,  1.907882,  1.167080,  0.025429},
				{ 16.639956,  5.636819,  0.108493, 47.299709,  0.379984,  0.000000}}},
	{ F,	-1,	{{  0.457649,  3.841561,  1.432771,  0.801876,  3.395041,  0.069525},
				{  0.917243,  5.507803,  0.164955, 51.076206, 15.821679,  0.000000}}},
	{ Na,	+1,	{{  3.148690,  4.073989,  0.767888,  0.995612,  0.968249,  0.045300},
				{  2.594987,  6.046925,  0.070139, 14.122657,  0.217037,  0.000000}}},
	{ Mg,	+2,	{{  3.062918,  4.135106,  0.853742,  1.036792,  0.852520,  0.058851},
				{  2.015803,  4.417941,  0.065307,  9.669710,  0.187818,  0.000000}}},
	{ Al,	+3,	{{  4.132015,  0.912049,  1.102425,  0.614876,  3.219136,  0.019397},
				{  3.528641,  7.378344,  0.133708,  0.039065,  1.644728,  0.000000}}},
	{ Si, atom_type_traits::kWKSFVal,
				{{  2.879033,  3.072960,  1.515981,  1.390030,  4.995051,  0.146030},
				{  1.239713, 38.706276,  0.081481, 93.616333,  2.770293,  0.000000}}},
	{ Si,	+4,	{{  3.676722,  3.828496,  1.258033,  0.419024,  0.720421,  0.097266},
				{  1.446851,  3.013144,  0.064397,  0.206254,  5.970222,  0.000000}}},
	{ Cl,	-1,	{{  1.061802,  7.139886,  6.524271,  2.355626, 35.829403,-34.916603},
				{  0.144727,  1.171795, 19.467655, 60.320301,  0.000436,  0.000000}}},
	{ K,	+1,	{{-17.609339,  1.494873,  7.150305, 10.899569, 15.808228,  0.257164},
				{ 18.840979,  0.053453,  0.812940, 22.264105, 14.351593,  0.000000}}},
	{ Ca,	+2,	{{  8.501441, 12.880483,  9.765095,  7.156669,  0.711160,-21.013187},
				{ 10.525848, -0.004033,  0.010692,  0.684443, 27.231771,  0.000000}}},
	{ Sc,	+3,	{{  7.104348,  1.511488,-53.669773, 38.404816, 24.532240,  0.118642},
				{  0.601957,  0.033386, 12.572138, 10.859736, 14.125230,  0.000000}}},
	{ Ti,	+2,	{{  7.040119,  1.496285,  9.657304,  0.006534,  1.649561,  0.150362},
				{  0.537072,  0.031914,  8.009958,201.800293, 24.039482,  0.000000}}},
	{ Ti,	+3,	{{ 36.587933,  7.230255, -9.086077,  2.084594, 17.294008,-35.111282},
				{  0.000681,  0.522262,  5.262317, 15.881716,  6.149805,  0.000000}}},
	{ Ti,	+4,	{{ 45.355537,  7.092900,  7.483858,-43.498817,  1.678915, -0.110628},
				{  9.252186,  0.523046, 13.082852, 10.193876,  0.023064,  0.000000}}},
	{ V,	+2,	{{  7.754356,  2.064100,  2.576998,  2.011404,  7.126177, -0.533379},
				{  7.066315,  0.014993,  7.066308, 22.055786,  0.467568,  0.000000}}},
	{ V,	+3,	{{  9.958480,  1.596350,  1.483442,-10.846044, 17.332867,  0.474921},
				{  6.763041,  0.056895, 17.750029,  0.328826,  0.388013,  0.000000}}},
	{ V,	+5,	{{ 15.575018,  8.448095,  1.612040, -9.721855,  1.534029,  0.552676},
				{  0.682708,  5.566640, 10.527077,  0.907961,  0.066667,  0.000000}}},
	{ Cr,	+2,	{{ 10.598877,  1.565858,  2.728280,  0.098064,  6.959321,  0.049870},
				{  6.151846,  0.023519, 17.432816, 54.002388,  0.426301,  0.000000}}},
	{ Cr,	+3,	{{  7.989310,  1.765079,  2.627125,  1.829380,  6.980908, -0.192123},
				{  6.068867,  0.018342,  6.068887, 16.309284,  0.420864,  0.000000}}},
	{ Mn,	+2,	{{ 11.287712, 26.042414,  3.058096,  0.090258,  7.088306,-24.566132},
				{  5.506225,  0.000774, 16.158575, 54.766354,  0.375580,  0.000000}}},
	{ Mn,	+3,	{{  6.926972,  2.081342, 11.128379,  2.375107, -0.419287, -0.093713},
				{  0.378315,  0.015054,  5.379957, 14.429586,  0.004939,  0.000000}}},
	{ Mn,	+4,	{{ 12.409131,  7.466993,  1.809947,-12.138477, 10.780248,  0.672146},
				{  0.300400,  0.112814, 12.520756,  0.168653,  5.173237,  0.000000}}},
	{ Fe,	+2,	{{ 11.776765, 11.165097,  3.533495,  0.165345,  7.036932, -9.676919},
				{  4.912232,  0.001748, 14.166556, 42.381958,  0.341324,  0.000000}}},
	{ Fe,	+3,	{{  9.721638, 63.403847,  2.141347,  2.629274,  7.033846,-61.930725},
				{  4.869297,  0.000293,  4.867602, 13.539076,  0.338520,  0.000000}}},
	{ Co,	+2,	{{  6.993840, 26.285812, 12.254289,  0.246114,  4.017407,-24.796852},
				{  0.310779,  0.000684,  4.400528, 35.741447, 12.536393,  0.000000}}},
	{ Co,	+3,	{{  6.861739,  2.678570, 12.281889,  3.501741, -0.179384, -1.147345},
				{  0.309794,  0.008142,  4.331703, 11.914167, 11.914167,  0.000000}}},
	{ Ni,	+2,	{{ 12.519017, 37.832058,  4.387257,  0.661552,  6.949072,-36.344471},
				{  3.933053,  0.000442, 10.449184, 23.860998,  0.283723,  0.000000}}},
	{ Ni,	+3,	{{ 13.579366,  1.902844, 12.859268,  3.811005, -6.838595, -0.317618},
				{  0.313140,  0.012621,  3.906407, 10.894311,  0.344379,  0.000000}}},
	{ Cu,	+1,	{{ 12.960763, 16.342150,  1.110102,  5.520682,  6.915452,-14.849320},
				{  3.576010,  0.000975, 29.523218, 10.114283,  0.261326,  0.000000}}},
	{ Cu,	+2,	{{ 11.895569, 16.344978,  5.799817,  1.048804,  6.789088,-14.878383},
				{  3.378519,  0.000924,  8.133653, 20.526524,  0.254741,  0.000000}}},
	{ Zn,	+2,	{{ 13.340772, 10.428857,  5.544489,  0.762295,  6.869172, -8.945248},
				{  3.215913,  0.001413,  8.542680, 21.891756,  0.239215,  0.000000}}},
	{ Ga,	+3,	{{ 13.123875, 35.288189,  6.126979,  0.611551,  6.724807,-33.875122},
				{  2.809960,  0.000323,  6.831534, 16.784311,  0.212002,  0.000000}}},
	{ Ge,	+4,	{{  6.876636,  6.779091,  9.969591,  3.135857,  0.152389,  1.086542},
				{  2.025174,  0.176650,  3.573822,  7.685848, 16.677574,  0.000000}}},
	{ Br,	-1,	{{ 17.714310,  6.466926,  6.947385,  4.402674, -0.697279,  1.152674},
				{  2.122554, 19.050768,  0.152708, 58.690361, 58.690372,  0.000000}}},
	{ Rb,	+1,	{{ 17.684320,  7.761588,  6.680874,  2.668883,  0.070974,  1.133263},
				{  1.710209, 14.919863,  0.128542, 31.654478,  0.128543,  0.000000}}},
	{ Sr,	+2,	{{ 17.694973,  1.275762,  6.154252,  9.234786,  0.515995,  1.125309},
				{  1.550888, 30.133041,  0.118774, 13.821799,  0.118774,  0.000000}}},
	{ Y,	+3,	{{ 46.660366, 10.369686,  4.623042,-62.170834, 17.471146, 19.023842},
				{ -0.019971, 13.180257,  0.176398, -0.016727,  1.467348,  0.000000}}},
	{ Zr,	+4,	{{  6.802956, 17.699253, 10.650647, -0.248108,  0.250338,  0.827902},
				{  0.096228,  1.296127, 11.240715, -0.219259, -0.219021,  0.000000}}},
	{ Nb,	+3,	{{ 17.714323,  1.675213,  7.483963,  8.322464, 11.143573, -8.339573},
				{  1.172419, 30.102791,  0.080255, -0.002983, 10.456687,  0.000000}}},
	{ Nb,	+5,	{{ 17.580206,  7.633277, 10.793497,  0.180884, 67.837921,-68.024780},
				{  1.165852,  0.078558,  9.507652, 31.621656, -0.000438,  0.000000}}},
	{ Mo,	+3,	{{  7.447050, 17.778122, 11.886068,  1.997905,  1.789626, -1.898764},
				{  0.072000,  1.073145,  9.834720, 28.221746, -0.011674,  0.000000}}},
	{ Mo,	+5,	{{  7.929879, 17.667669, 11.515987,  0.500402, 77.444084,-78.056595},
				{  0.068856,  1.068064,  9.046229, 26.558945, -0.000473,  0.000000}}},
	{ Mo,	+6,	{{ 34.757683,  9.653037,  6.584769,-18.628115,  2.490594,  1.141916},
				{  1.301770,  7.123843,  0.094097,  1.617443, 12.335434,  0.000000}}},
	{ Ru,	+3,	{{ 17.894758, 13.579529, 10.729251,  2.474095, 48.227997,-51.905243},
				{  0.902827,  8.740579,  0.045125, 24.764954, -0.001699,  0.000000}}},
	{ Ru,	+4,	{{ 17.845776, 13.455084, 10.229087,  1.653524, 14.059795,-17.241762},
				{  0.901070,  8.482392,  0.045972, 23.015272, -0.004889,  0.000000}}},
	{ Rh,	+3,	{{ 17.758621, 14.569813,  5.298320,  2.533579,  0.879753,  0.960843},
				{  0.841779,  8.319533,  0.069050, 23.709131,  0.069050,  0.000000}}},
	{ Rh,	+4,	{{ 17.716188, 14.446654,  5.185801,  1.703448,  0.989992,  0.959941},
				{  0.840572,  8.100647,  0.068995, 22.357307,  0.068995,  0.000000}}},
	{ Pd,	+2,	{{  6.122282, 15.651012,  3.513508,  9.060790,  8.771199,  0.879336},
				{  0.062424,  8.018296, 24.784275,  0.776457,  0.776457,  0.000000}}},
	{ Pd,	+4,	{{  6.152421,-96.069023, 31.622141, 81.578255, 17.801403,  0.915874},
				{  0.063951, 11.090354, 13.466152,  9.758302,  0.783014,  0.000000}}},
	{ Ag,	+1,	{{  6.091192,  4.019526, 16.948174,  4.258638, 13.889437,  0.785127},
				{  0.056305,  0.719340,  7.758938, 27.368349,  0.719340,  0.000000}}},
	{ Ag,	+2,	{{  6.401808, 48.699802,  4.799859,-32.332523, 16.356710,  1.068247},
				{  0.068167,  0.942270, 20.639496,  1.100365,  6.883131,  0.000000}}},
	{ Cd,	+2,	{{  6.093711, 43.909691, 17.041306,-39.675117, 17.958918,  0.664795},
				{  0.050624,  8.654143, 15.621396, 11.082067,  0.667591,  0.000000}}},
	{ In,	+3,	{{  6.206277, 18.497746,  3.078131, 10.524613,  7.401234,  0.293677},
				{  0.041357,  6.605563, 18.792250,  0.608082,  0.608082,  0.000000}}},
	{ Sn,	+2,	{{  6.353672,  4.770377, 14.672025,  4.235959, 18.002131, -0.042519},
				{  0.034720,  6.167891,  6.167879, 29.006456,  0.561774,  0.000000}}},
	{ Sn,	+4,	{{ 15.445732,  6.420892,  4.562980,  1.713385, 18.033537, -0.172219},
				{  6.280898,  0.033144,  6.280899, 17.983601,  0.557980,  0.000000}}},
	{ Sb,	+3,	{{ 10.189171, 57.461918, 19.356573,  4.862206,-45.394096,  1.516108},
				{  0.089485,  0.375256,  5.357987, 22.153736,  0.297768,  0.000000}}},
	{ Sb,	+5,	{{ 17.920622,  6.647932, 12.724075,  1.555545,  7.600591, -0.445371},
				{  0.522315,  0.029487,  5.718210, 16.433775,  5.718204,  0.000000}}},
	{ I,	-1,	{{ 20.010330, 17.835524,  8.104130,  2.231118,  9.158548, -3.341004},
				{  4.565931,  0.444266, 32.430672, 95.149040,  0.014906,  0.000000}}},
	{ Cs,	+1,	{{ 19.939056, 24.967621, 10.375884,  0.454243, 17.660248,-19.394306},
				{  3.770511,  0.004040, 25.311275, 76.537766,  0.384730,  0.000000}}},
	{ Ba,	+2,	{{ 19.750200, 17.513683, 10.884892,  0.321585, 65.149834,-59.618172},
				{  3.430748,  0.361590, 21.358307, 70.309402,  0.001418,  0.000000}}},
	{ La,	+3,	{{ 19.688887, 17.345703, 11.356296,  0.099418, 82.358124,-76.846909},
				{  3.146211,  0.339586, 18.753832, 90.345459,  0.001072,  0.000000}}},
	{ Ce,	+3,	{{ 26.593231, 85.866432, -6.677695, 12.111847, 17.401903,-80.313423},
				{  3.280381,  0.001012,  4.313575, 17.868504,  0.326962,  0.000000}}},
	{ Ce,	+4,	{{ 17.457533, 25.659941, 11.691037, 19.695251,-16.994749, -3.515096},
				{  0.311812, -0.003793, 16.568687,  2.886395, -0.008931,  0.000000}}},
	{ Pr,	+3,	{{ 20.879841, 36.035797, 12.135341,  0.283103, 17.167803,-30.500784},
				{  2.870897,  0.002364, 16.615236, 53.909359,  0.306993,  0.000000}}},
	{ Pr,	+4,	{{ 17.496082, 21.538509, 20.403114, 12.062211, -7.492043, -9.016722},
				{  0.294457, -0.002742,  2.772886, 15.804613, -0.013556,  0.000000}}},
	{ Nd,	+3,	{{ 17.120077, 56.038139, 21.468307, 10.000671,  2.905866,-50.541992},
				{  0.291295,  0.001421,  2.743681, 14.581367, 22.485098,  0.000000}}},
	{ Pm,	+3,	{{ 22.221066, 17.068142, 12.805423,  0.435687, 52.238770,-46.767181},
				{  2.635767,  0.277039, 14.927315, 45.768017,  0.001455,  0.000000}}},
	{ Sm,	+3,	{{ 15.618565, 19.538092, 13.398946, -4.358811, 24.490461, -9.714854},
				{  0.006001,  0.306379, 14.979594,  0.748825,  2.454492,  0.000000}}},
	{ Eu,	+2,	{{ 23.899035, 31.657497, 12.955752,  1.700576, 16.992199,-26.204315},
				{  2.467332,  0.002230, 13.625002, 35.089481,  0.253136,  0.000000}}},
	{ Eu,	+3,	{{ 17.758327, 33.498665, 24.067188, 13.436883, -9.019134,-19.768026},
				{  0.244474, -0.003901,  2.487526, 14.568011, -0.015628,  0.000000}}},
	{ Gd,	+3,	{{ 24.344999, 16.945311, 13.866931,  0.481674, 93.506378,-88.147179},
				{  2.333971,  0.239215, 12.982995, 43.876347,  0.000673,  0.000000}}},
	{ Tb,	+3,	{{ 24.878252, 16.856016, 13.663937,  1.279671, 39.271294,-33.950317},
				{  2.223301,  0.227290, 11.812528, 29.910065,  0.001527,  0.000000}}},
	{ Dy,	+3,	{{ 16.864344, 90.383461, 13.675473,  1.687078, 25.540651,-85.150650},
				{  0.216275,  0.000593, 11.121207, 26.250975,  2.135930,  0.000000}}},
	{ Ho,	+3,	{{ 16.837524, 63.221336, 13.703766,  2.061602, 26.202621,-58.026505},
				{  0.206873,  0.000796, 10.500283, 24.031883,  2.055060,  0.000000}}},
	{ Er,	+3,	{{ 16.810127, 22.681061, 13.864114,  2.294506, 26.864477,-17.513460},
				{  0.198293,  0.002126,  9.973341, 22.836388,  1.979442,  0.000000}}},
	{ Tm,	+3,	{{ 16.787500, 15.350905, 14.182357,  2.299111, 27.573771,-10.192087},
				{  0.190852,  0.003036,  9.602934, 22.526880,  1.912862,  0.000000}}},
	{ Yb,	+2,	{{ 28.443794, 16.849527, 14.165081,  3.445311, 28.308853,-23.214935},
				{  1.863896,  0.183811,  9.225469, 23.691355,  0.001463,  0.000000}}},
	{ Yb,	+3,	{{ 28.191629, 16.828087, 14.167848,  2.744962, 23.171774,-18.103676},
				{  1.842889,  0.182788,  9.045957, 20.799847,  0.001759,  0.000000}}},
	{ Lu,	+3,	{{ 28.828693, 16.823227, 14.247617,  3.079559, 25.647667,-20.626528},
				{  1.776641,  0.175560,  8.575531, 19.693701,  0.001453,  0.000000}}},
	{ Hf,	+4,	{{ 29.267378, 16.792543, 14.785310,  2.184128, 23.791996,-18.820383},
				{  1.697911,  0.168313,  8.190025, 18.277578,  0.001431,  0.000000}}},
	{ Ta,	+5,	{{ 29.539469, 16.741854, 15.182070,  1.642916, 16.437447,-11.542459},
				{  1.612934,  0.160460,  7.654408, 17.070732,  0.001858,  0.000000}}},
	{ W,	+6,	{{ 29.729357, 17.247808, 15.184488,  1.154652,  0.739335,  3.945157},
				{  1.501648,  0.140803,  6.880573, 14.299601, 14.299618,  0.000000}}},
	{ Os,	+4,	{{ 17.113485, 15.792370, 23.342392,  4.090271,  7.671292,  3.988390},
				{  0.131850,  7.288542,  1.389307, 19.629425,  1.389307,  0.000000}}},
	{ Ir,	+3,	{{ 31.537575, 16.363338, 15.597141,  5.051404,  1.436935,  4.009459},
				{  1.334144,  7.451918,  0.127514, 21.705648,  0.127515,  0.000000}}},
	{ Ir,	+4,	{{ 30.391249, 16.146996, 17.019068,  4.458904,  0.975372,  4.006865},
				{  1.328519,  7.181766,  0.127337, 19.060146,  1.328519,  0.000000}}},
	{ Pt,	+2,	{{ 31.986849, 17.249048, 15.269374,  5.760234,  1.694079,  4.032512},
				{  1.281143,  7.625512,  0.123571, 24.190826,  0.123571,  0.000000}}},
	{ Pt,	+4,	{{ 41.932713, 16.339224, 17.653894,  6.012420,-12.036877,  4.094551},
				{  1.111409,  6.466086,  0.128917, 16.954155,  0.778721,  0.000000}}},
	{ Au,	+1,	{{ 32.124306, 16.716476, 16.814100,  7.311565,  0.993064,  4.040792},
				{  1.216073,  7.165378,  0.118715, 20.442486, 53.095985,  0.000000}}},
	{ Au,	+3,	{{ 31.704271, 17.545767, 16.819551,  5.522640,  0.361725,  4.042679},
				{  1.215561,  7.220506,  0.118812, 20.050970,  1.215562,  0.000000}}},
	{ Hg,	+1,	{{ 28.866837, 19.277540, 16.776051,  6.281459,  3.710289,  4.068430},
				{  1.173967,  7.583842,  0.115351, 29.055994,  1.173968,  0.000000}}},
	{ Hg,	+2,	{{ 32.411079, 18.690371, 16.711773,  9.974835, -3.847611,  4.052869},
				{  1.162980,  7.329806,  0.114518, 22.009489, 22.009493,  0.000000}}},
	{ Tl,	+1,	{{ 32.295044, 16.570049, 17.991013,  1.535355,  7.554591,  4.054030},
				{  1.101544,  0.110020,  6.528559, 52.495068, 20.338634,  0.000000}}},
	{ Tl,	+3,	{{ 32.525639, 19.139185, 17.100321,  5.891115, 12.599463, -9.256075},
				{  1.094966,  6.900992,  0.103667, 18.489614, -0.001401,  0.000000}}},
	{ Pb,	+2,	{{ 27.392647, 16.496822, 19.984501,  6.813923,  5.233910,  4.065623},
				{  1.058874,  0.106305,  6.708123, 24.395554,  1.058874,  0.000000}}},
	{ Pb,	+4,	{{ 32.505657, 20.014240, 14.645661,  5.029499,  1.760138,  4.044678},
				{  1.047035,  6.670321,  0.105279, 16.525040,  0.105279,  0.000000}}},
	{ Bi,	+3,	{{ 32.461437, 19.438683, 16.302486,  7.322662,  0.431704,  4.043703},
				{  0.997930,  6.038867,  0.101338, 18.371586, 46.361046,  0.000000}}},
	{ Bi,	+5,	{{ 16.734028, 20.580494,  9.452623, 61.155834,-34.041023,  4.113663},
				{  0.105076,  4.773282, 11.762162,  1.211775,  1.619408,  0.000000}}},
	{ Ra,	+2,	{{  4.986228, 32.474945, 21.947443, 11.800013, 10.807292,  3.956572},
				{  0.082597,  0.791468,  4.608034, 24.792431,  0.082597,  0.000000}}},
	{ Ac,	+3,	{{ 15.584983, 32.022125, 21.456327,  0.757593, 12.341252,  3.838984},
				{  0.077438,  0.739963,  4.040735, 47.525002, 19.406845,  0.000000}}},
	{ Th,	+4,	{{ 15.515445, 32.090691, 13.996399, 12.918157,  7.635514,  3.831122},
				{  0.074499,  0.711663,  3.871044, 18.596891,  3.871044,  0.000000}}},
	{ U,	+3,	{{ 15.360309, 32.395657, 21.961290,  1.325894, 14.251453,  3.706622},
				{  0.067815,  0.654643,  3.643409, 39.604965, 16.330570,  0.000000}}},
	{ U,	+4,	{{ 15.355091, 32.235306,  0.557745, 14.396367, 21.751173,  3.705863},
				{  0.067789,  0.652613, 42.354237, 15.908239,  3.553231,  0.000000}}},
	{ U,	+6,	{{ 15.333844, 31.770849, 21.274414, 13.872636,  0.048519,  3.700591},
				{  0.067644,  0.646384,  3.317894, 14.650250, 75.339699,  0.000000}}},
	{ Np,	+3,	{{ 15.378152, 32.572132, 22.206125,  1.413295, 14.828381,  3.603370},
				{  0.064613,  0.631420,  3.561936, 37.875511, 15.546129,  0.000000}}},
	{ Np,	+4,	{{ 15.373926, 32.423019, 21.969994,  0.662078, 14.969350,  3.603039},
				{  0.064597,  0.629658,  3.476389, 39.438942, 15.135764,  0.000000}}},
	{ Np,	+6,	{{ 15.359986, 31.992825, 21.412458,  0.066574, 14.568174,  3.600942},
				{  0.064528,  0.624505,  3.253441, 67.658318, 13.980832,  0.000000}}},
	{ Pu,	+3,	{{ 15.356004, 32.769127, 22.680210,  1.351055, 15.416232,  3.428895},
				{  0.060590,  0.604663,  3.491509, 37.260635, 14.981921,  0.000000}}},
	{ Pu,	+4,	{{ 15.416219, 32.610569, 22.256662,  0.719495, 15.518152,  3.480408},
				{  0.061456,  0.607938,  3.411848, 37.628792, 14.464360,  0.000000}}},
	{ Pu,	+6,	{{ 15.436506, 32.289719, 14.726737, 15.012391,  7.024677,  3.502325},
				{  0.061815,  0.606541,  3.245363, 13.616438,  3.245364,  0.000000}}}
};

SFDataArrayElement kELSFData[] = {
	{H, 0,	{{	0.034900,	0.120100,	0.197000,	0.057300,	0.119500	},
			{	0.534700,	3.586700,	12.347100,	18.952499,	38.626900	}}},
	{He, 0,	{{	0.031700,	0.083800,	0.152600,	0.133400,	0.016400	},
			{	0.250700,	1.475100,	4.493800,	12.664600,	31.165300	}}},
	{Li, 0,	{{	0.075000,	0.224900,	0.554800,	1.495400,	0.935400	},
			{	0.386400,	2.938300,	15.382900,	53.554501,	138.733704	}}},
	{Be, 0,	{{	0.078000,	0.221000,	0.674000,	1.386700,	0.692500	},
			{	0.313100,	2.238100,	10.151700,	30.906099,	78.327301	}}},
	{B, 0,	{{	0.090900,	0.255100,	0.773800,	1.213600,	0.460600	},
			{	0.299500,	2.115500,	8.381600,	24.129200,	63.131401	}}},
	{C, 0,	{{	0.089300,	0.256300,	0.757000,	1.048700,	0.357500	},
			{	0.246500,	1.710000,	6.409400,	18.611300,	50.252300	}}},
	{N, 0,	{{	0.102200,	0.321900,	0.798200,	0.819700,	0.171500	},
			{	0.245100,	1.748100,	6.192500,	17.389400,	48.143101	}}},
	{O, 0,	{{	0.097400,	0.292100,	0.691000,	0.699000,	0.203900	},
			{	0.206700,	1.381500,	4.694300,	12.710500,	32.472599	}}},
	{F, 0,	{{	0.108300,	0.317500,	0.648700,	0.584600,	0.142100	},
			{	0.205700,	1.343900,	4.278800,	11.393200,	28.788099	}}},
	{Ne, 0,	{{	0.126900,	0.353500,	0.558200,	0.467400,	0.146000	},
			{	0.220000,	1.377900,	4.020300,	9.493400,	23.127800	}}},
	{Na, 0,	{{	0.214200,	0.685300,	0.769200,	1.658900,	1.448200	},
			{	0.333400,	2.344600,	10.083000,	48.303699,	138.270004	}}},
	{Mg, 0,	{{	0.231400,	0.686600,	0.967700,	2.188200,	1.133900	},
			{	0.327800,	2.272000,	10.924100,	39.289799,	101.974800	}}},
	{Al, 0,	{{	0.239000,	0.657300,	1.201100,	2.558600,	1.231200	},
			{	0.313800,	2.106300,	10.416300,	34.455200,	98.534401	}}},
	{Si, 0,	{{	0.251900,	0.637200,	1.379500,	2.508200,	1.050000	},
			{	0.307500,	2.017400,	9.674600,	29.374399,	80.473198	}}},
	{P, 0,	{{	0.254800,	0.610600,	1.454100,	2.320400,	0.847700	},
			{	0.290800,	1.874000,	8.517600,	24.343399,	63.299599	}}},
	{S, 0,	{{	0.249700,	0.562800,	1.389900,	2.186500,	0.771500	},
			{	0.268100,	1.671100,	7.026700,	19.537701,	50.388802	}}},
	{Cl, 0,	{{	0.244300,	0.539700,	1.391900,	2.019700,	0.662100	},
			{	0.246800,	1.524200,	6.153700,	16.668699,	42.308601	}}},
	{Ar, 0,	{{	0.238500,	0.501700,	1.342800,	1.889900,	0.607900	},
			{	0.228900,	1.369400,	5.256100,	14.092800,	35.536098	}}},
	{K, 0,	{{	0.411500,	1.403100,	2.278400,	2.674200,	2.216200	},
			{	0.370300,	3.387400,	13.102900,	68.959198,	194.432907	}}},
	{Ca, 0,	{{	0.405400,	1.388000,	2.160200,	3.753200,	2.206300	},
			{	0.349900,	3.099100,	11.960800,	53.935299,	142.389206	}}},
	{Sc, 0,	{{	0.378700,	1.218100,	2.059400,	3.261800,	2.387000	},
			{	0.313300,	2.585600,	9.581300,	41.768799,	116.728203	}}},
	{Ti, 0,	{{	0.382500,	1.259800,	2.000800,	3.061700,	2.069400	},
			{	0.304000,	2.486300,	9.278300,	39.075100,	109.458298	}}},
	{V, 0,	{{	0.387600,	1.275000,	1.910900,	2.831400,	1.897900	},
			{	0.296700,	2.378000,	8.798100,	35.952801,	101.720100	}}},
	{Cr, 0,	{{	0.404600,	1.369600,	1.894100,	2.080000,	1.219600	},
			{	0.298600,	2.395800,	9.140600,	37.470100,	113.712097	}}},
	{Mn, 0,	{{	0.379600,	1.209400,	1.781500,	2.542000,	1.593700	},
			{	0.269900,	2.045500,	7.472600,	31.060400,	91.562202	}}},
	{Fe, 0,	{{	0.394600,	1.272500,	1.703100,	2.314000,	1.479500	},
			{	0.271700,	2.044300,	7.600700,	29.971399,	86.226501	}}},
	{Co, 0,	{{	0.411800,	1.316100,	1.649300,	2.193000,	1.283000	},
			{	0.274200,	2.037200,	7.720500,	29.968000,	84.938301	}}},
	{Ni, 0,	{{	0.386000,	1.176500,	1.545100,	2.073000,	1.381400	},
			{	0.247800,	1.766000,	6.310700,	25.220400,	74.314598	}}},
	{Cu, 0,	{{	0.431400,	1.320800,	1.523600,	1.467100,	0.856200	},
			{	0.269400,	1.922300,	7.347400,	28.989201,	90.624603	}}},
	{Zn, 0,	{{	0.428800,	1.264600,	1.447200,	1.829400,	1.093400	},
			{	0.259300,	1.799800,	6.750000,	25.586000,	73.528397	}}},
	{Ga, 0,	{{	0.481800,	1.403200,	1.656100,	2.460500,	1.105400	},
			{	0.282500,	1.978500,	8.754600,	32.523800,	98.552299	}}},
	{Ge, 0,	{{	0.465500,	1.301400,	1.608800,	2.699800,	1.300300	},
			{	0.264700,	1.792600,	7.607100,	26.554100,	77.523804	}}},
	{As, 0,	{{	0.451700,	1.222900,	1.585200,	2.795800,	1.263800	},
			{	0.249300,	1.643600,	6.815400,	22.368099,	62.039001	}}},
	{Se, 0,	{{	0.447700,	1.167800,	1.584300,	2.808700,	1.195600	},
			{	0.240500,	1.544200,	6.323100,	19.461000,	52.023300	}}},
	{Br, 0,	{{	0.479800,	1.194800,	1.869500,	2.695300,	0.820300	},
			{	0.250400,	1.596300,	6.965300,	19.849199,	50.323299	}}},
	{Kr, 0,	{{	0.454600,	1.099300,	1.769600,	2.706800,	0.867200	},
			{	0.230900,	1.427900,	5.944900,	16.675200,	42.224300	}}},
	{Rb, 0,	{{	1.016000,	2.852800,	3.546600,	-7.780400,	12.114800	},
			{	0.485300,	5.092500,	25.785101,	130.451508,	138.677505	}}},
	{Sr, 0,	{{	0.670300,	1.492600,	3.336800,	4.460000,	3.150100	},
			{	0.319000,	2.228700,	10.350400,	52.329102,	151.221603	}}},
	{Y, 0,	{{	0.689400,	1.547400,	3.245000,	4.212600,	2.976400	},
			{	0.318900,	2.290400,	10.006200,	44.077099,	125.012001	}}},
	{Zr, 0,	{{	0.671900,	1.468400,	3.166800,	3.955700,	2.892000	},
			{	0.303600,	2.124900,	8.923600,	36.845798,	108.204903	}}},
	{Nb, 0,	{{	0.612300,	1.267700,	3.034800,	3.384100,	2.368300	},
			{	0.270900,	1.768300,	7.248900,	27.946501,	98.562401	}}},
	{Mo, 0,	{{	0.677300,	1.479800,	3.178800,	3.082400,	1.838400	},
			{	0.292000,	2.060600,	8.112900,	30.533600,	100.065804	}}},
	{Tc, 0,	{{	0.708200,	1.639200,	3.199300,	3.432700,	1.871100	},
			{	0.297600,	2.210600,	8.524600,	33.145599,	96.637703	}}},
	{Ru, 0,	{{	0.673500,	1.493400,	3.096600,	2.725400,	1.559700	},
			{	0.277300,	1.971600,	7.324900,	26.689100,	90.558098	}}},
	{Rh, 0,	{{	0.641300,	1.369000,	2.985400,	2.695200,	1.543300	},
			{	0.258000,	1.772100,	6.385400,	23.254900,	85.151703	}}},
	{Pd, 0,	{{	0.590400,	1.177500,	2.651900,	2.287500,	0.868900	},
			{	0.232400,	1.501900,	5.159100,	15.542800,	46.821301	}}},
	{Ag, 0,	{{	0.637700,	1.379000,	2.829400,	2.363100,	1.455300	},
			{	0.246600,	1.697400,	5.765600,	20.094299,	76.737198	}}},
	{Cd, 0,	{{	0.636400,	1.424700,	2.780200,	2.597300,	1.788600	},
			{	0.240700,	1.682300,	5.658800,	20.721901,	69.110901	}}},
	{In, 0,	{{	0.676800,	1.658900,	2.774000,	3.183500,	2.132600	},
			{	0.252200,	1.854500,	6.293600,	25.145700,	84.544800	}}},
	{Sn, 0,	{{	0.722400,	1.961000,	2.716100,	3.560300,	1.897200	},
			{	0.265100,	2.060400,	7.301100,	27.549299,	81.334900	}}},
	{Sb, 0,	{{	0.710600,	1.924700,	2.614900,	3.832200,	1.889900	},
			{	0.256200,	1.964600,	6.885200,	24.764799,	68.916801	}}},
	{Te, 0,	{{	0.694700,	1.869000,	2.535600,	4.001300,	1.895500	},
			{	0.245900,	1.854200,	6.441100,	22.173000,	59.220600	}}},
	{I, 0,	{{	0.704700,	1.948400,	2.594000,	4.152600,	1.505700	},
			{	0.245500,	1.863800,	6.763900,	21.800699,	56.439499	}}},
	{Xe, 0,	{{	0.673700,	1.790800,	2.412900,	4.210000,	1.705800	},
			{	0.230500,	1.689000,	5.821800,	18.392799,	47.249599	}}},
	{Cs, 0,	{{	1.270400,	3.801800,	5.661800,	0.920500,	4.810500	},
			{	0.435600,	4.205800,	23.434200,	136.778305,	171.756104	}}},
	{Ba, 0,	{{	0.904900,	2.607600,	4.849800,	5.160300,	4.738800	},
			{	0.306600,	2.436300,	12.182100,	54.613499,	161.997803	}}},
	{La, 0,	{{	0.840500,	2.386300,	4.613900,	5.151400,	4.794900	},
			{	0.279100,	2.141000,	10.340000,	41.914799,	132.020401	}}},
	{Ce, 0,	{{	0.855100,	2.391500,	4.577200,	5.027800,	4.511800	},
			{	0.280500,	2.120000,	10.180800,	42.063301,	130.989304	}}},
	{Pr, 0,	{{	0.909600,	2.531300,	4.526600,	4.637600,	4.369000	},
			{	0.293900,	2.247100,	10.826600,	48.884201,	147.602005	}}},
	{Nd, 0,	{{	0.880700,	2.418300,	4.444800,	4.685800,	4.172500	},
			{	0.280200,	2.083600,	10.035700,	47.450600,	146.997604	}}},
	{Pm, 0,	{{	0.947100,	2.546300,	4.352300,	4.478900,	3.908000	},
			{	0.297700,	2.227600,	10.576200,	49.361900,	145.358002	}}},
	{Sm, 0,	{{	0.969900,	2.583700,	4.277800,	4.457500,	3.598500	},
			{	0.300300,	2.244700,	10.648700,	50.799400,	146.417892	}}},
	{Eu, 0,	{{	0.869400,	2.241300,	3.919600,	3.969400,	4.549800	},
			{	0.265300,	1.859000,	8.399800,	36.739700,	125.708900	}}},
	{Gd, 0,	{{	0.967300,	2.470200,	4.114800,	4.497200,	3.209900	},
			{	0.290900,	2.101400,	9.706700,	43.426998,	125.947403	}}},
	{Tb, 0,	{{	0.932500,	2.367300,	3.879100,	3.967400,	3.799600	},
			{	0.276100,	1.951100,	8.929600,	41.593700,	131.012207	}}},
	{Dy, 0,	{{	0.950500,	2.370500,	3.821800,	4.047100,	3.445100	},
			{	0.277300,	1.946900,	8.886200,	43.093800,	133.139603	}}},
	{Ho, 0,	{{	0.924800,	2.242800,	3.618200,	3.791000,	3.791200	},
			{	0.266000,	1.818300,	7.965500,	33.112900,	101.813904	}}},
	{Er, 0,	{{	1.037300,	2.482400,	3.655800,	3.892500,	3.005600	},
			{	0.294400,	2.079700,	9.415600,	45.805599,	132.772003	}}},
	{Tm, 0,	{{	1.007500,	2.378700,	3.544000,	3.693200,	3.175900	},
			{	0.281600,	1.948600,	8.716200,	41.841999,	125.031998	}}},
	{Yb, 0,	{{	1.034700,	2.391100,	3.461900,	3.655600,	3.005200	},
			{	0.285500,	1.967900,	8.761900,	42.330399,	125.649902	}}},
	{Lu, 0,	{{	0.992700,	2.243600,	3.355400,	3.781300,	3.099400	},
			{	0.270100,	1.807300,	7.811200,	34.484901,	103.352600	}}},
	{Hf, 0,	{{	1.029500,	2.291100,	3.411000,	3.949700,	2.492500	},
			{	0.276100,	1.862500,	8.096100,	34.271198,	98.529503	}}},
	{Ta, 0,	{{	1.019000,	2.229100,	3.409700,	3.925200,	2.267900	},
			{	0.269400,	1.796200,	7.694400,	31.094200,	91.108902	}}},
	{W, 0,	{{	0.985300,	2.116700,	3.357000,	3.798100,	2.279800	},
			{	0.256900,	1.674500,	7.009800,	26.923401,	81.390999	}}},
	{Re, 0,	{{	0.991400,	2.085800,	3.453100,	3.881200,	1.852600	},
			{	0.254800,	1.651800,	6.884500,	26.723400,	81.721497	}}},
	{Os, 0,	{{	0.981300,	2.032200,	3.366500,	3.623500,	1.974100	},
			{	0.248700,	1.597300,	6.473700,	23.281700,	70.925400	}}},
	{Ir, 0,	{{	1.019400,	2.064500,	3.442500,	3.491400,	1.697600	},
			{	0.255400,	1.647500,	6.596600,	23.226900,	70.027199	}}},
	{Pt, 0,	{{	0.914800,	1.809600,	3.213400,	3.295300,	1.575400	},
			{	0.226300,	1.381300,	5.324300,	17.598700,	60.017101	}}},
	{Au, 0,	{{	0.967400,	1.891600,	3.399300,	3.052400,	1.260700	},
			{	0.235800,	1.471200,	5.675800,	18.711901,	61.528599	}}},
	{Hg, 0,	{{	1.003300,	1.946900,	3.439600,	3.154800,	1.418000	},
			{	0.241300,	1.529800,	5.800900,	19.452000,	60.575298	}}},
	{Tl, 0,	{{	1.068900,	2.103800,	3.603900,	3.492700,	1.828300	},
			{	0.254000,	1.671500,	6.350900,	23.153099,	78.709900	}}},
	{Pb, 0,	{{	1.089100,	2.186700,	3.616000,	3.803100,	1.899400	},
			{	0.255200,	1.717400,	6.513100,	23.917000,	74.703903	}}},
	{Bi, 0,	{{	1.100700,	2.230600,	3.568900,	4.154900,	2.038200	},
			{	0.254600,	1.735100,	6.494800,	23.646400,	70.377998	}}},
	{Po, 0,	{{	1.156800,	2.435300,	3.645900,	4.406400,	1.717900	},
			{	0.264800,	1.878600,	7.174900,	25.176600,	69.282097	}}},
	{At, 0,	{{	1.090900,	2.197600,	3.383100,	4.670000,	2.127700	},
			{	0.246600,	1.670700,	6.019700,	20.765699,	57.266300	}}},
	{Rn, 0,	{{	1.075600,	2.163000,	3.317800,	4.885200,	2.048900	},
			{	0.240200,	1.616900,	5.764400,	19.456800,	52.500900	}}},
	{Fr, 0,	{{	1.428200,	3.508100,	5.676700,	4.196400,	3.894600	},
			{	0.318300,	2.688900,	13.481600,	54.386600,	200.832108	}}},
	{Ra, 0,	{{	1.312700,	3.124300,	5.298800,	5.389100,	5.413300	},
			{	0.288700,	2.289700,	10.827600,	43.538898,	145.610901	}}},
	{Ac, 0,	{{	1.312800,	3.102100,	5.338500,	5.961100,	4.756200	},
			{	0.286100,	2.250900,	10.528700,	41.779598,	128.297302	}}},
	{Th, 0,	{{	1.255300,	2.917800,	5.086200,	6.120600,	4.712200	},
			{	0.270100,	2.063600,	9.305100,	34.597698,	107.919998	}}},
	{Pa, 0,	{{	1.321800,	3.144400,	5.437100,	5.644400,	4.010700	},
			{	0.282700,	2.225000,	10.245400,	41.116199,	124.444901	}}},
	{U, 0,	{{	1.338200,	3.204300,	5.455800,	5.483900,	3.634200	},
			{	0.283800,	2.245200,	10.251900,	41.725101,	124.902298	}}},
	{Np, 0,	{{	1.519300,	4.005300,	6.532700,	-0.140200,	6.748900	},
			{	0.321300,	2.820600,	14.887800,	68.910301,	81.725700	}}},
	{Pu, 0,	{{	1.351700,	3.293700,	5.321300,	4.646600,	3.571400	},
			{	0.281300,	2.241800,	9.995200,	42.793900,	132.173904	}}},
	{Am, 0,	{{	1.213500,	2.796200,	4.754500,	4.573100,	4.478600	},
			{	0.248300,	1.843700,	7.542100,	29.384100,	112.457901	}}},
	{Cm, 0,	{{	1.293700,	3.110000,	5.039300,	4.754600,	3.503100	},
			{	0.263800,	2.034100,	8.710100,	35.299198,	109.497200	}}},
	{Bk, 0,	{{	1.291500,	3.102300,	4.930900,	4.600900,	3.466100	},
			{	0.261100,	2.002300,	8.437700,	34.155899,	105.891098	}}},
	{Cf, 0,	{{	1.208900,	2.739100,	4.348200,	4.004700,	4.649700	},
			{	0.242100,	1.748700,	6.726200,	23.215300,	80.310799	}}}
};

} // namespace data

// --------------------------------------------------------------------
// atom_type_traits

atom_type_traits::atom_type_traits(const std::string& symbol)
	: m_info(nullptr)
{
	for (auto& i: data::kKnownAtoms)
	{
		if (cif::iequals(i.symbol, symbol))
		{
			m_info = &i;
			break;
		}
	}

	if (symbol == "X")
		m_info = &data::kKnownAtoms[0];
	
	if (m_info == nullptr)
		throw std::invalid_argument("Not a known element: " + symbol);
}

atom_type_traits::atom_type_traits(atom_type t)
{
	if (t < H or t >= data::kKnownAtomsCount)
		throw std::invalid_argument("atomType out of range");
	
	m_info = &data::kKnownAtoms[t];
	
	assert(m_info->type == t);
}

bool atom_type_traits::is_element(const std::string& symbol)
{
	bool result = false;
	
	for (auto& i: data::kKnownAtoms)
	{
		if (cif::iequals(i.symbol, symbol))
		{
			result = true;
			break;
		}
	}
	
	return result;
}

bool atom_type_traits::is_metal(const std::string& symbol)
{
	bool result = false;
	
	for (auto& i: data::kKnownAtoms)
	{
		if (cif::iequals(i.symbol, symbol))
		{
			result = i.metal;
			break;
		}
	}
	
	return result;
}

bool atom_type_traits::has_sf(int charge) const
{
	auto type = m_info->type;
	if (type == D)
		type = H;

	bool result = false;

	for (auto& sf: data::kWKSFData)
	{
		if (sf.symbol == type and sf.charge == charge)
		{
			result = true;
			break;
		}
	}

	return result;
}

auto atom_type_traits::wksf(int charge) const -> const SFData&
{
	auto type = m_info->type;
	if (type == D)
		type = H;

	for (auto& sf: data::kWKSFData)
	{
		if (sf.symbol == type and sf.charge == charge)
			return sf.sf;
	}

	if (charge != 0)
	{
		// Oops, not found. Fall back to zero charge and see if we can use that

		if (cif::VERBOSE > 0)
			std::cerr << "No scattering factor found for " << name() << " with charge " << charge << " will try to fall back to zero charge...\n";

		for (auto& sf: data::kWKSFData)
		{
			if (sf.symbol == type and sf.charge == 0)
				return sf.sf;
		}
	}

	throw std::out_of_range("No scattering factor found for " + name() + std::to_string(charge));
}
	
auto atom_type_traits::elsf() const -> const SFData&
{
	auto type = m_info->type;
	if (type == D)
		type = H;

	for (auto& sf: data::kELSFData)
	{
		if (sf.symbol == type)
			return sf.sf;
	}

	throw std::invalid_argument("No scattering factor found for " + name());
}

// ionic radii

float atom_type_traits::crystal_ionic_radius(int charge) const
{
	float result = kNA;

	if (charge >= -3 and charge <= 8)
	{
		for (auto &r : data::kCrystalIonicRadii)
		{
			if (r.type != m_info->type)
				continue;
			
			result = r.radii[charge < 0 ? charge + 3 : charge + 2] / 100.0f;
			break;
		}
	}

	return result;
}

float atom_type_traits::effective_ionic_radius(int charge) const
{
	float result = kNA;

	if (charge >= -3 and charge <= 8)
	{
		for (auto &r : data::kEffectiveIonicRadii)
		{
			if (r.type != m_info->type)
				continue;
			
			result = r.radii[charge < 0 ? charge + 3 : charge + 2] / 100.0f;
			break;
		}
	}

	return result;
}

}
