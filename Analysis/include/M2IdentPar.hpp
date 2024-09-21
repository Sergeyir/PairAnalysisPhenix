// $HEADER$
//------------------------------------------------------------------------------------------------
//                            M2IdentPar variables declaration
//------------------------------------------------------------------------------------------------
// M2IdentPar - m2 identification parameters
//
// ** Code for use in PHENIX related projects **
//
// Author: Sergei Antsupov
// Email: antsupov0124@gmail.com
//
/**
 * Header for storing global variables for m2 identification of pions, kaons, and protons
 **/
//------------------------------------------------------------------------------------------------

#ifndef M2_IDENT_PAR_HPP
#define M2_IDENT_PAR_HPP

// Info on parameters indices of sigma functions
// mean = pol1
// sigma:
// [0] = sigma_{alpha}
// [1] = sigma_{ms}
// [2] = sigma_{t}
// [3] = K1
// [4] = L
//
// Means were fitted with 1st order polynomial

#ifdef RUN7AUAU

static const double M2_MEAN_PAR_PION_TOFE[7] = {0.0174898, 0.00427226};
static const double M2_MEAN_PAR_KAON_TOFE[7] = {0.245808, 1.15555e-05};
static const double M2_MEAN_PAR_PROTON_TOFE[7] = {0.889027, -0.0127422};
static const double M2_MEAN_PAR_APION_TOFE[7] = {0.0175598, -0.00436738};
static const double M2_MEAN_PAR_AKAON_TOFE[7] = {0.245897, -0.000129887};
static const double M2_MEAN_PAR_APROTON_TOFE[7] = {0.894373, 0.0147789};
static const double M2_SIGMA_PAR_TOFE[5] = {0.835, 1., 120, 104, 5.2};

static const double M2_MEAN_PAR_TOFW_PION[7] = {0.0195054, 0.00083337};
static const double M2_MEAN_PAR_TOFW_KAON[7] = {0.246532, -0.00120144};
static const double M2_MEAN_PAR_TOFW_PROTON[7] = {0.889054, -0.00509407};
static const double M2_MEAN_PAR_TOFW_APION[7] = {0.0194966, -0.000710241};
static const double M2_MEAN_PAR_TOFW_AKAON[7] = {0.246375, 0.00113778};
static const double M2_MEAN_PAR_TOFW_APROTON[7] = {0.887917, 0.00325808};
static const double M2_SIGMA_PAR_TOFW[5] = {0.9, 1., 75.9, 104, 4.95};

//[sector][parameter]
static const double M2_MEAN_PAR_EMCALE_PION[2][2] = 
	{{0.0193227, 0.00266133}, {0.0197781, 0.00250303}};
static const double M2_MEAN_PAR_EMCALE_KAON[2][2] = 
	{{0.271721, -0.027387}, {0.283546, -0.0420572}};
static const double M2_MEAN_PAR_EMCALE_PROTON[2][2] = 
	{{0.866507, 0.0149589}, {0.883922, 0.0113606}};
static const double M2_MEAN_PAR_EMCALE_APION[2][2] = 
	{{0.0170196, -0.0064267}, {0.0174645, -0.00567766}};
static const double M2_MEAN_PAR_EMCALE_AKAON[2][2] = 
	{{0.249519, 0.00512605}, {0.252869, 0.00774439}};
static const double M2_MEAN_PAR_EMCALE_APROTON[2][2] = 
	{{0.875881, 0.028631}, {0.887292, 0.0373331}};
static const double M2_SIGMA_PAR_EMCALE[2][5] = 
	{{0.835, 2.3, 500., 104., 5.22}, {0.835, 2.4, 520., 104., 5.22}};

static const double M2_MEAN_PAR_EMCALW_PION[4][2] = 
	{{0.02156, -0.000862084}, {0.0212173, 3.7905e-05}, 
	{0.0212057, 0.00102381}, {0.0209914, 0.00158534}};
static const double M2_MEAN_PAR_EMCALW_KAON[4][2] = 
	{{0.277186, -0.0260212}, {0.281834, -0.0293181}, 
	{0.288341, -0.0306334}, {0.286166, -0.025915}};
static const double M2_MEAN_PAR_EMCALW_PROTON[4][2] = 
	{{0.902183, 0.00659663}, {0.917261, -0.00396265}, 
	{0.942533, -0.00469523}, {0.934547, 0.00766581}};
static const double M2_MEAN_PAR_EMCALW_APION[4][2] = 
	{{0.0191554, -0.00321943}, {0.0188638, -0.00438432}, 
	{0.0188607, -0.00465492}, {0.0189317, -0.00487683}};
static const double M2_MEAN_PAR_EMCALW_AKAON[4][2] = 
	{{0.255676, 0.00350508}, {0.262825, 0.00987294}, 
	{0.266624, 0.00799149}, {0.266556, 0.00601535}};
static const double M2_MEAN_PAR_EMCALW_APROTON[4][2] = 
	{{0.907618, 0.0351048}, {0.932191, 0.0537718}, 
	{0.960642, 0.0565557}, {0.952011, 0.0432118}};
static const double M2_SIGMA_PAR_EMCALW[4][5] = 
	{{0.9, 1.8, 480., 104., 5.22}, {0.9, 1.7, 480., 104., 5.22}, 
	{0.9, 1.5, 480., 104., 5.22}, {0.9, 2., 500., 104., 5.22}};

#endif /* RUN7AUAU */

#endif /* M2_IDENT_PAR_HPP */
