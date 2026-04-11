/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 FieldLine
//
///////////////////////////////////////////////////////////////////////////////

#include "FieldLine.h"

#pragma warning(disable : 4251 4100 4244 4101)

//////////////////////////////////////////////////////////////////////////
// definition of class FieldLine
//////////////////////////////////////////////////////////////////////////
vtCFieldLine::vtCFieldLine(CVectorField* pField):
m_nNumSeeds(0),
m_integrationOrder(FOURTH),
m_timeDir(FORWARD),
m_fInitTime((double)0.0),
m_fStepTime((double)0.0),
m_fDurationTime((double)0.0),
m_fLowerAngleAccuracy((double)0.99),
m_fUpperAngleAccuracy((double)0.999),
m_fStationaryCutoff((double)0.00001),
m_nMaxsize(MAX_LENGTH),
m_nSaveInterval(1),
m_lastFieldLineStepCount(0),
m_fInitStepSize(1.0),
m_pField(pField)
{
}

vtCFieldLine::~vtCFieldLine(void)
{
	releaseSeedMemory();
}

//////////////////////////////////////////////////////////////////////////
// release the memory allocated to seeds
//////////////////////////////////////////////////////////////////////////
void vtCFieldLine::releaseSeedMemory(void)
{
	vtListParticleIter pIter = m_lSeeds.begin();
	for( ; pIter != m_lSeeds.end(); ++pIter )
	{
		vtParticleInfo* thisPart = *pIter;
		delete thisPart;
	}
	m_lSeeds.erase(m_lSeeds.begin(), m_lSeeds.end() );
	m_lSeedIds.erase(m_lSeedIds.begin(), m_lSeedIds.end() );
	m_nNumSeeds = 0;
}

//////////////////////////////////////////////////////////////////////////
// Integrate along a field line using the 2nd order Euler-Cauchy 
// predictor-corrector method. This routine is used for both steady and
// unsteady vector fields. 
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::euler_cauchy(TIME_DIR, TIME_DEP,double*, double)
{
	return 1;
}

//////////////////////////////////////////////////////////////////////////
// Integrate along a field line using the 2th order Runge-Kutta method.
// This routine is used for both steady and unsteady vector fields. 
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::runge_kutta2(TIME_DIR time_dir, TIME_DEP time_dep, 
							   PointInfo& ci, 
							   double* t, double dt,
							   int* cachedLowT)
{
	int istat = 0;

	return istat;
}

//////////////////////////////////////////////////////////////////////////
// Integrate along a field line using the 4th order Runge-Kutta method.
// This routine is used for both steady and unsteady vector fields. 
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::runge_kutta4(TIME_DIR time_dir, 
							   TIME_DEP time_dep, 
							   PointInfo& ci, 
							   double* t,			// initial time
							   double dt,			// stepsize
							   int* cachedLowT)
{
	int i, istat;
	VECTOR3 pt0;
	VECTOR3 vel;
	VECTOR3 k1, k2, k3;
	VECTOR3 pt;
	int fromCell;

	pt = ci.phyCoord;
	// 1st step of the Runge-Kutta scheme
	istat = m_pField->at_phys(ci.fromCell, pt, ci, *t, vel, cachedLowT);
	if ( istat != 1 )
		return OUT_OF_BOUND;

	for( i=0; i<3; i++ )
	{
		pt0[i] = pt[i];
		k1[i] = time_dir*dt*vel[i];
		pt[i] = pt0[i]+k1[i]*(double)0.5;
	}

	// 2nd step of the Runge-Kutta scheme
	fromCell = ci.inCell;
	if ( time_dep  == UNSTEADY)
		*t += (double)0.5*time_dir*dt;
	
	istat=m_pField->at_phys(fromCell, pt, ci, *t, vel, cachedLowT);
	if ( istat!= 1 )
	{
		ci.phyCoord = pt;
		return OUT_OF_BOUND;
	}
	
	for( i=0; i<3; i++ )
	{
		k2[i] = time_dir*dt*vel[i];
		pt[i] = pt0[i]+k2[i]*(double)0.5;
	}

	// 3rd step of the Runge-Kutta scheme
	fromCell = ci.inCell;
	istat=m_pField->at_phys(fromCell, pt, ci, *t, vel, cachedLowT);
	if ( istat != 1 )
	{
		ci.phyCoord = pt;
		return OUT_OF_BOUND;
	}
	
	for( i=0; i<3; i++ )
	{
		k3[i] = time_dir*dt*vel[i];
		pt[i] = pt0[i]+k3[i];
	}

	//    4th step of the Runge-Kutta scheme
	if ( time_dep  == UNSTEADY)
		*t += (double)0.5*time_dir*dt;
	
	fromCell = ci.inCell;
	istat=m_pField->at_phys(fromCell, pt, ci, *t, vel, cachedLowT);
	if ( istat != 1 )
	{
		ci.phyCoord = pt;
		return OUT_OF_BOUND;
	}
	
	for( i=0; i<3; i++ )
	{
		pt[i] = pt0[i]+(k1[i]+(double)2.0*(k2[i]+k3[i])+time_dir*dt*vel[i])/(double)6.0;
	}
	ci.phyCoord = pt;

	return( istat );
}

//////////////////////////////////////////////////////////////////////////
// Geometric helper: one Euler step on the MPAS-O sphere.
// Takes horizontal velocity h_vel and vertical velocity v_vel already
// evaluated at pt_src, and advances the particle by time fdt.
// Returns false only if the particle goes below the planetary centre.
// When h_vel is zero the particle is left in place (returns true).
//////////////////////////////////////////////////////////////////////////
bool vtCFieldLine::geodesic_step(TIME_DIR   time_dir,
								  VECTOR3 pt_src, double r_src,
								  VECTOR3 h_vel,  double v_vel,
								  double fdt,
								  VECTOR3& pt_dst, double& r_dst)
{
	double vel_mag = h_vel.GetMag() * fdt;
	if (vel_mag == 0.0) {
		pt_dst = pt_src;
		r_dst  = r_src;
		return true;
	}
	double r_new = r_src + time_dir * v_vel * fdt;
	if (r_new <= 0.0) return false;
	double  omega    = vel_mag / r_src;
	VECTOR3 normal   = time_dir * cross(pt_src, h_vel);
	MATRIX3 rotate_m = rotate_matrix_axis(normal, omega);
	pt_dst = rotate_m * pt_src;
	pt_dst.Normalize();
	pt_dst.scale(r_new);
	r_dst = r_new;
	return true;
}

//////////////////////////////////////////////////////////////////////////
// Integrate along a field line using the Euler Runge-Kutta method.
// This routine is used for both steady and unsteady vector fields. 
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::MPASO_euler(TIME_DIR time_dir, 
							   TIME_DEP time_dep, 
							   PointInfo& ci, 
							   double* t,			// initial time
							   double dt,			// stepsize
							   int* cachedLowT)
{
	int istat;
	VECTOR4 vel;
	VECTOR3 horizontal_vec;
	double vertical_vec;
	VECTOR3 pt;
	int fromCell;

	pt = ci.phyCoord;
	// Euler
	fromCell = ci.inCell;
	istat = m_pField->at_phys(fromCell, pt, ci, *t, vel, cachedLowT);
	if ( istat != 1 )
		return OUT_OF_BOUND;

	horizontal_vec[0] = vel[0];
	horizontal_vec[1] = vel[1];
	horizontal_vec[2] = vel[2];
	vertical_vec = vel[3];
	double radius = pt.GetMag();
	if (radius <= 0.0) return OUT_OF_BOUND;  // degenerate position

	VECTOR3 pt_dst;
	double  r_dst;
	if (!geodesic_step(time_dir, pt, radius, horizontal_vec, vertical_vec, dt, pt_dst, r_dst))
		return OUT_OF_BOUND;

	ci.phyCoord = pt_dst;
	ci.fromCell = ci.inCell;
	if (time_dep == UNSTEADY)
		*t += time_dir * dt;
	return istat;
}

//////////////////////////////////////////////////////////////////////////
// 4th-order Runge-Kutta integration on the MPAS-O sphere.
// Mirrors MPASO_euler but evaluates velocity at four stages and combines
// them with the classical RK4 weights (1/6, 1/3, 1/3, 1/6).
// Horizontal motion uses geodesic rotation; vertical uses radial shift.
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::MPASO_rk4(TIME_DIR time_dir,
							   TIME_DEP time_dep,
							   PointInfo& ci,
							   double* t,
							   double dt,
							   int* cachedLowT)
{
	const double speedup = 100.0;
	int istat;
	VECTOR4 vel;

	VECTOR3 pt0 = ci.phyCoord;
	double  r0  = pt0.GetMag();
	if (r0 <= 0.0) return OUT_OF_BOUND;

	// --- Stage 1: evaluate k1 at (pt0, t) and cap dt ---
	PointInfo ci_tmp = ci;
	istat = m_pField->at_phys(ci.inCell, pt0, ci_tmp, *t, vel, cachedLowT);
	if (istat != 1) return OUT_OF_BOUND;

	VECTOR3 k1_h(vel[0], vel[1], vel[2]);
	double  k1_v = vel[3];
	int     cell0 = ci_tmp.inCell;

	// --- Stage 2: k2 at (pt0 + k1 * dt/2, t + dt/2) ---
	VECTOR3 pt1; double r1;
	if (!geodesic_step(time_dir, pt0, r0, k1_h, k1_v, 0.5 * dt, pt1, r1)) return OUT_OF_BOUND;
	double t1 = (time_dep == UNSTEADY) ? *t + 0.5 * time_dir * dt : *t;
	ci_tmp = ci; ci_tmp.phyCoord = pt1;
	istat = m_pField->at_phys(cell0, pt1, ci_tmp, t1, vel, cachedLowT);
	if (istat != 1) { ci.phyCoord = pt1; return OUT_OF_BOUND; }
	VECTOR3 k2_h(vel[0], vel[1], vel[2]);
	double  k2_v = vel[3];

	// --- Stage 3: k3 at (pt0 + k2 * dt/2, t + dt/2) ---
	VECTOR3 pt2; double r2;
	if (!geodesic_step(time_dir, pt0, r0, k2_h, k2_v, 0.5 * dt, pt2, r2)) return OUT_OF_BOUND;
	ci_tmp = ci; ci_tmp.phyCoord = pt2;
	istat = m_pField->at_phys(cell0, pt2, ci_tmp, t1, vel, cachedLowT);
	if (istat != 1) { ci.phyCoord = pt2; return OUT_OF_BOUND; }
	VECTOR3 k3_h(vel[0], vel[1], vel[2]);
	double  k3_v = vel[3];

	// --- Stage 4: k4 at (pt0 + k3 * dt, t + dt) ---
	VECTOR3 pt3; double r3;
	if (!geodesic_step(time_dir, pt0, r0, k3_h, k3_v, dt, pt3, r3)) return OUT_OF_BOUND;
	double t2 = (time_dep == UNSTEADY) ? *t + time_dir * dt : *t;
	ci_tmp = ci; ci_tmp.phyCoord = pt3;
	istat = m_pField->at_phys(cell0, pt3, ci_tmp, t2, vel, cachedLowT);
	if (istat != 1) { ci.phyCoord = pt3; return OUT_OF_BOUND; }
	VECTOR3 k4_h(vel[0], vel[1], vel[2]);
	double  k4_v = vel[3];

	// --- RK4 weighted average: (k1 + 2*k2 + 2*k3 + k4) / 6 ---
	VECTOR3 v_avg_h;
	for (int i = 0; i < 3; i++)
		v_avg_h[i] = (k1_h[i] + 2.0*k2_h[i] + 2.0*k3_h[i] + k4_h[i]) / 6.0;
	double v_avg_v = (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v) / 6.0;

	// --- Apply one rotation step with the averaged velocity ---
	VECTOR3 pt_final; double r_final;
	if (!geodesic_step(time_dir, pt0, r0, v_avg_h, v_avg_v, dt, pt_final, r_final)) return OUT_OF_BOUND;

	ci.phyCoord = pt_final;
	ci.fromCell = ci.inCell;

	if (time_dep == UNSTEADY)
		*t += time_dir * dt;

	return istat;
}

//////////////////////////////////////////////////////////////////////////
// aptive step size
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::adapt_step(const VECTOR3& p2,
							 const VECTOR3& p1,
							 const VECTOR3& p0,
							 double dt_estimate,
							 double* dt)
{
	double angle;
	VECTOR3 p2p1 = p2 - p1;
	VECTOR3 p1p0 = p1 - p0;
		
	angle = (double)acos(dot(p1p0, p2p1)/(p1p0.GetMag() * p2p1.GetMag()))*(double)RAD_TO_DEG;
	if(angle > m_fUpperAngleAccuracy)
		*dt = (*dt) * (double)0.5;
	else
		if(angle < m_fLowerAngleAccuracy)
		{
			*dt = (*dt) * (double)2.0;
			if(*dt >= m_fMaxStepSize)
				*dt = m_fMaxStepSize;
		}
	
	return true;
}

//////////////////////////////////////////////////////////////////////////
// initialize seeds
//////////////////////////////////////////////////////////////////////////
void vtCFieldLine::setSeedPoints(VECTOR3* points, int numPoints, double t,
		int64_t *seedIds)
{
	int i, res;
	VECTOR3 nodeData;

	if(points == NULL)
		return;
	
	// if the rake size has changed, forget the previous seed points
	if( m_nNumSeeds != numPoints )
		releaseSeedMemory();
	
	if ( seedIds != NULL )
	{
		for (i = 0 ; i < numPoints ; i++ )
		{
			m_lSeedIds.push_back( seedIds[i] );
		}
	}
	
	if( m_nNumSeeds == 0 )
	{
		for(i = 0 ; i < numPoints ; i++ )
		{
			vtParticleInfo* newParticle = new vtParticleInfo;
			newParticle->m_pointInfo.phyCoord = points[i];
			newParticle->m_fStartTime = t;
			newParticle->m_fCurrentTime = t;
			newParticle->m_cachedLowT = -1;
			newParticle->ptId = i;

			// query the field in order to get the starting
			// cell interpolant for the seed point
			// which was passed to us without any prior information
			res = m_pField->at_phys(-1, points[i], newParticle->m_pointInfo, t, nodeData, &newParticle->m_cachedLowT);
// 			newParticle->itsValidFlag =  (res == 1) ? 1 : 0 ;
			newParticle->itsValidFlag =  1;
			m_lSeeds.push_back( newParticle );
		}
	} 
	else
	{
		vtListParticleIter sIter = m_lSeeds.begin();

		for(i = 0 ; sIter != m_lSeeds.end() ; ++sIter, ++i )
		{
			vtParticleInfo* thisSeed = *sIter;

			// set the new location for this seed point
			thisSeed->m_pointInfo.phyCoord = points[i];
			thisSeed->m_fStartTime = t;
			thisSeed->m_fCurrentTime = t;
			thisSeed->m_cachedLowT = -1;
			res = m_pField->at_phys(-1, points[i], thisSeed->m_pointInfo, t, nodeData, &thisSeed->m_cachedLowT);
			thisSeed->itsValidFlag =  (res == 1) ? 1 : 0 ;
		}    
	}
	m_nNumSeeds = numPoints;
}




//////////////////////////////////////////////////////////////////////////
// initialize seeds with possibly different start times in t array 
//////////////////////////////////////////////////////////////////////////
void vtCFieldLine::setSeedPoints(VECTOR3* points, int numPoints, 
				 double* t)
{
	int i, res;
	VECTOR3 nodeData;

	if(points == NULL)
		return;

	// if the rake size has changed, forget the previous seed points
	if( m_nNumSeeds != numPoints )
		releaseSeedMemory();

	if( m_nNumSeeds == 0 )
	{
		for(i = 0 ; i < numPoints ; i++ )
		{
			vtParticleInfo* newParticle = new vtParticleInfo;

			newParticle->m_pointInfo.phyCoord = points[i];
			newParticle->m_fStartTime = t[i];
			newParticle->m_fCurrentTime = t[i];
			newParticle->m_cachedLowT = -1;
			newParticle->ptId = i;

			// query the field in order to get the starting
			// cell interpolant for the seed point
			// which was passed to us without any prior information

			res = m_pField->at_phys(-1, points[i], newParticle->m_pointInfo, t[i], nodeData, &newParticle->m_cachedLowT);
// 			newParticle->itsValidFlag =  (res == 1) ? 1 : 0 ;
			newParticle->itsValidFlag =  1;
			m_lSeeds.push_back( newParticle );
		}
	} 
	else
	{
		vtListParticleIter sIter = m_lSeeds.begin();

		for(i = 0 ; sIter != m_lSeeds.end() ; ++sIter, ++i )
		{
			vtParticleInfo* thisSeed = *sIter;

			// set the new location for this seed point
			thisSeed->m_pointInfo.phyCoord = points[i];
			thisSeed->m_fStartTime = t[i];
			thisSeed->m_fCurrentTime = t[i];
			thisSeed->m_cachedLowT = -1;
			res = m_pField->at_phys(-1, points[i], thisSeed->m_pointInfo, t[i], nodeData, &thisSeed->m_cachedLowT);
			thisSeed->itsValidFlag =  (res == 1) ? 1 : 0 ;
		}    
	}

	m_nNumSeeds = numPoints;
}



//////////////////////////////////////////////////////////////////////////
// initialize seeds with possibly different start times 
//////////////////////////////////////////////////////////////////////////
void vtCFieldLine::setSeedPoints(VECTOR4* points, int numPoints, 
		int64_t *seedIds)
{
	int i, res;
	VECTOR3 nodeData;

	if(points == NULL)
		return;

	// if the rake size has changed, forget the previous seed points
	if( m_nNumSeeds != numPoints )
		releaseSeedMemory();

	if ( seedIds != NULL )
	{
		for (i = 0 ; i < numPoints ; i++ )
		{
			m_lSeedIds.push_back( seedIds[i] );
		}
	}

	if( m_nNumSeeds == 0 )
	{
		for(i = 0 ; i < numPoints ; i++ )
		{
			vtParticleInfo* newParticle = new vtParticleInfo;

			VECTOR3 pos(points[i][0], points[i][1], points[i][2]);
			newParticle->m_pointInfo.phyCoord = pos;
			newParticle->m_fStartTime = points[i][3];
			newParticle->m_fCurrentTime = points[i][3];
			newParticle->m_cachedLowT = -1;
			newParticle->ptId = i;

			// query the field in order to get the starting
			// cell interpolant for the seed point
			// which was passed to us without any prior information

			res = m_pField->at_phys(-1, pos, newParticle->m_pointInfo, points[i][3], nodeData, &newParticle->m_cachedLowT);
// 			newParticle->itsValidFlag =  (res == 1) ? 1 : 0 ;
			newParticle->itsValidFlag =  1;
			m_lSeeds.push_back( newParticle );
		}
	} 
	else
	{
		vtListParticleIter sIter = m_lSeeds.begin();

		for(i = 0 ; sIter != m_lSeeds.end() ; ++sIter, ++i )
		{
			vtParticleInfo* thisSeed = *sIter;

			VECTOR3 pos(points[i][0], points[i][1], points[i][2]); 
			// set the new location for this seed point
			thisSeed->m_pointInfo.phyCoord = pos;
			thisSeed->m_fStartTime = points[i][3];
			thisSeed->m_fCurrentTime = points[i][3];
			thisSeed->m_cachedLowT = -1;
			res = m_pField->at_phys(-1, pos, thisSeed->m_pointInfo, points[i][3], nodeData, &thisSeed->m_cachedLowT);
			thisSeed->itsValidFlag =  (res == 1) ? 1 : 0 ;
		}    
	}

	m_nNumSeeds = numPoints;
}
