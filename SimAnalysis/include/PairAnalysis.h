#pragma once

//struct for storing id of particles
struct id_container
{
	int tof;
	int emc;
	int pc2;
	int pc3;
	int orig_id;
};

bool IsTOF2PID(id_container *Id1, id_container *Id2)
{
	if (Id1->tof == Id1->orig_id && Id2->tof == Id2->orig_id) return true;
	return false;
}

bool IsEMC2PID(id_container *Id1, id_container *Id2)
{
	if (Id1->emc == Id1->orig_id && Id2->emc == Id2->orig_id) return true;

	return false;
}

bool IsEMCnoPID(id_container *Id1, id_container *Id2)
{
	if (Id1->emc != PartId.junk && Id2->emc != PartId.junk) return true;

	return false;
}

bool Is1TOFandEMC2PID(id_container *Id1, id_container *Id2)
{
	if 
	(
		IsTOF2PID(Id1, Id2) ||
		(Id1->emc == Id1->orig_id && Id2->tof == Id2->orig_id) || 
		(Id1->tof == Id1->orig_id && Id2->emc == Id2->orig_id)
	) return true;

	return false;
}

bool Is2PID(id_container *Id1, id_container *Id2)
{
	if (IsEMC2PID(Id1, Id2) || Is1TOFandEMC2PID(Id1, Id2)) return true;
	
	return false;
}

bool Is1PID(id_container *Id1, id_container *Id2)
{
	if 
	(
		(Id1->tof == Id1->orig_id && 
		(
			Id2->tof != PartId.junk || 
			Id2->emc != PartId.junk || 
			Id2->pc2 != PartId.junk ||
			Id2->pc3 != PartId.junk)
		) || 
		(Id2->tof == Id2->orig_id && 
		(
			Id1->tof != PartId.junk || 
			Id1->emc != PartId.junk || 
			Id1->pc2 != PartId.junk ||
			Id1->pc3 != PartId.junk)
		) 
	) return true;
	
	return false;
}

bool IsnoPID(id_container *Id1, id_container *Id2)
{
	return true;
}

double GetPairPt(const double *pp1, const double *pp2)
{
	const double pt = sqrt((pp1[0]+pp2[0])*(pp1[0]+pp2[0]) + (pp1[1]+pp2[1])*(pp1[1]+pp2[1]));
	return pt;
}

double GetMass(const double *pp1, const double *pp2, const double m1, const double m2)
{
	const double pm1 = pp1[0]*pp1[0] + pp1[1]*pp1[1] + pp1[2]*pp1[2];
	const double pm2 = pp2[0]*pp2[0] + pp2[1]*pp2[1] + pp2[2]*pp2[2];

	if( pm1 <= 0. || pm2 <= 0. ) return -999;

	const double e1 = sqrt(pow(m1, 2) + pm1);
	const double e2 = sqrt(pow(m2, 2) + pm2);
	const double es = e1 + e2;

	double ps[3];
	
	ps[0] = pp1[0] + pp2[0];
	ps[1] = pp1[1] + pp2[1];
	ps[2] = pp1[2] + pp2[2];
	
	const double mass = sqrt(es*es - ps[0]*ps[0] - ps[1]*ps[1] - ps[2]*ps[2]);
	
	return mass;
}

bool IsOneArmCut(const double phi_pip, const double phi_pim)
{
	if ( !((phi_pip<1.5&&phi_pim<1.5) || (phi_pip>1.5&&phi_pim>1.5) )) return true;
	return false;
}

bool IsGhostCut(double dzed, double dphi, double dalpha)
{
	//pc1
	if (fabs(dzed) < 6.0 && fabs(dphi - (0.13*dalpha)) < 0.015) return true;
	//x1x2_1
	if (fabs(dphi - (0.04*dalpha)) < 0.015) return true;
	//x1x2_2
	if (fabs(dphi - (-0.065*dalpha)) < 0.015) return true;
	return false;
}

bool noCut()
{
	return true;
}
