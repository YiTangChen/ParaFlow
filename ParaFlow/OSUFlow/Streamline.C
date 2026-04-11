/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vector Field
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 Streamlines
//
///////////////////////////////////////////////////////////////////////////////

#include "FieldLine.h"

#pragma warning(disable : 4251 4100 4244 4101)

FILE* fDebugOut;
//////////////////////////////////////////////////////////////////////////
// definition of class FieldLine
//////////////////////////////////////////////////////////////////////////
vtCStreamLine::vtCStreamLine(CVectorField* pField):
vtCFieldLine(pField),
m_itsTraceDir(BACKWARD_AND_FORWARD),
m_fPsuedoTime(0.0),
m_fCurrentTime(0.0)
{
#ifdef DEBUG
	fDebugOut = fopen(".\\DebugInfo\\debug.txt", "w");
#endif
}

vtCStreamLine::~vtCStreamLine(void)
{
#ifdef DEBUG
	fclose(fDebugOut);
#endif
}

//////////////////////////////////////////////////////////////////////////
// Set the advection direction
//////////////////////////////////////////////////////////////////////////
void vtCStreamLine::setForwardTracing(int enabled)
{
	switch( m_itsTraceDir )
	{
	case BACKWARD_DIR:
		if(enabled)
			m_itsTraceDir = BACKWARD_AND_FORWARD;
		break;
	case FORWARD_DIR:
		if(!enabled)
			m_itsTraceDir = OFF;
		break;
	case BACKWARD_AND_FORWARD:
		if(!enabled)
			m_itsTraceDir = BACKWARD_DIR;
		break;
	case OFF:
		if (enabled)
			m_itsTraceDir = FORWARD_DIR;
		break;
	}
}

void vtCStreamLine::setBackwardTracing(int enabled)
{
	switch( m_itsTraceDir )
	{
	  case BACKWARD_DIR:
		  if(!enabled)
			  m_itsTraceDir = OFF;
		  break;
	  case FORWARD_DIR:
		  if(enabled)
			  m_itsTraceDir = BACKWARD_AND_FORWARD;
		  break;
	  case BACKWARD_AND_FORWARD:
		  if(!enabled)
			  m_itsTraceDir = FORWARD_DIR;
		  break;
	  case OFF:
		  if (enabled)
			  m_itsTraceDir = BACKWARD_DIR;
		  break;
	}
}

int vtCStreamLine::getForwardTracing(void)
{
	if( m_itsTraceDir==FORWARD_DIR || m_itsTraceDir==BACKWARD_AND_FORWARD )
		return 1;
	else
		return 0;
}

int vtCStreamLine::getBackwardTracing(void)
{
	if( m_itsTraceDir==BACKWARD_DIR || m_itsTraceDir==BACKWARD_AND_FORWARD )
		return 1;
	else
		return 0;
}

//////////////////////////////////////////////////////////////////////////
// Compute streamlines.
// output
//	listSeedTraces: For each seed, return a list keeping the trace it
//					advects
//	listSeedIds: For each seed, a unique id may be associated with it. 
//					These unique ids will populate the list (defaults to NULL)
//////////////////////////////////////////////////////////////////////////
void vtCStreamLine::execute(const void* userData, 
			    list<vtListSeedTrace*>& listSeedTraces,
					list<int64_t> *listSeedIds)
{
	m_fCurrentTime = *(double *)userData;
	computeStreamLine(userData, listSeedTraces, listSeedIds);
}

void vtCStreamLine::execute(const void* userData, 
				list<vtListSeedTrace*>& listSeedTraces,
				vector<int>& fromCells,
				vector<int>& toCells,
				MPASOReader* mpaso_reader,
				list<int64_t> *listSeedIds)
{
	m_fCurrentTime = *(double *)userData;
	computeStreamLine(userData, listSeedTraces, fromCells, toCells, mpaso_reader, listSeedIds);
}

void vtCStreamLine::computeStreamLine(const void* userData,
				      list<vtListSeedTrace*>& listSeedTraces,
							list<int64_t> *listSeedIds)
{

	vtListParticleIter sIter;
	list<int64_t>::iterator sIdIter;
	int iterCnt = 0;

	#if	0	// MOD-BY-LEETEN 02/07/2011-FROM:
	for(sIter = m_lSeeds.begin(), sIdIter = m_lSeedIds.begin(); 
			sIter != m_lSeeds.end(); ++sIter, ++sIdIter)
	#else	// MOD-BY-LEETEN 02/07/2011-TO:
	if( !m_lSeedIds.empty() ) 
			sIdIter = m_lSeedIds.begin();
	for(sIter = m_lSeeds.begin(); sIter != m_lSeeds.end(); ++sIter)
	#endif	// MOD-BY-LEETEN 02/07/2011-END
	{
		vtParticleInfo* thisSeed = *sIter;

		if(thisSeed->itsValidFlag == 1)			// valid seed
		{

			if(m_itsTraceDir & BACKWARD_DIR)
			{
				vtListSeedTrace* backTrace;
				backTrace = new vtListSeedTrace;
				computeFieldLine(BACKWARD,m_integrationOrder, STEADY, *backTrace, thisSeed->m_pointInfo);
				listSeedTraces.push_back(backTrace);
				if (listSeedIds != NULL)
					(*listSeedIds).push_back(*sIdIter);
			}
			if(m_itsTraceDir & FORWARD_DIR)
			{
				vtListSeedTrace* forwardTrace;
				forwardTrace = new vtListSeedTrace;
				computeFieldLine(FORWARD,m_integrationOrder, STEADY, *forwardTrace, thisSeed->m_pointInfo);
				listSeedTraces.push_back(forwardTrace);
				if (listSeedIds != NULL)
					(*listSeedIds).push_back(*sIdIter);
			}
			iterCnt++;
		}
		// ADD-BY-LEETEN 02/07/2011-BEGIN
		if( !m_lSeedIds.empty() ) 
			sIdIter++;
		// ADD-BY-LEETEN 02/07/2011-END
	}
}

void vtCStreamLine::computeStreamLine(const void* userData,
				      list<vtListSeedTrace*>& listSeedTraces,
					  vector<int>& fromCells,
					  vector<int>& toCells,
					  MPASOReader* mpaso_reader,
					  list<int64_t> *listSeedIds)
{

	vtListParticleIter sIter;
	list<int64_t>::iterator sIdIter;
	int iterCnt = 0;

	m_traceActualStepCounts.clear();

	#if	0	// MOD-BY-LEETEN 02/07/2011-FROM:
	for(sIter = m_lSeeds.begin(), sIdIter = m_lSeedIds.begin();
			sIter != m_lSeeds.end(); ++sIter, ++sIdIter)
	#else	// MOD-BY-LEETEN 02/07/2011-TO:
	if( !m_lSeedIds.empty() )
			sIdIter = m_lSeedIds.begin();
	for(sIter = m_lSeeds.begin(); sIter != m_lSeeds.end(); ++sIter)
	#endif	// MOD-BY-LEETEN 02/07/2011-END
	{
		vtParticleInfo* thisSeed = *sIter;

		if(thisSeed->itsValidFlag == 1)			// valid seed
		{
			int fromCell = -1;
			if(fromCells.size()>iterCnt)
				fromCell = fromCells[iterCnt];
			if(m_itsTraceDir & BACKWARD_DIR)
			{
				vtListSeedTrace* backTrace;
				backTrace = new vtListSeedTrace;
				computeFieldLine(BACKWARD,m_integrationOrder, STEADY, *backTrace, thisSeed->m_pointInfo, fromCell, toCells);
				listSeedTraces.push_back(backTrace);
				m_traceActualStepCounts.push_back(m_lastFieldLineStepCount);
				if (listSeedIds != NULL)
					(*listSeedIds).push_back(*sIdIter);
			}
			if(m_itsTraceDir & FORWARD_DIR)
			{
				vtListSeedTrace* forwardTrace;
				forwardTrace = new vtListSeedTrace;
				computeFieldLine(FORWARD,m_integrationOrder, STEADY, *forwardTrace, thisSeed->m_pointInfo, fromCell, toCells);
				listSeedTraces.push_back(forwardTrace);
				m_traceActualStepCounts.push_back(m_lastFieldLineStepCount);
				if (listSeedIds != NULL)
					(*listSeedIds).push_back(*sIdIter);
			}
			iterCnt++;
		}
		// ADD-BY-LEETEN 02/07/2011-BEGIN
		if( !m_lSeedIds.empty() )
			sIdIter++;
		// ADD-BY-LEETEN 02/07/2011-END
	}
}

int vtCStreamLine::computeFieldLine( TIME_DIR time_dir,
									 INTEG_ORD integ_ord,
									 TIME_DEP time_dep, 
									 vtListSeedTrace& seedTrace,
									 PointInfo& seedInfo)
{
	int count = 0, istat;
	PointInfo thisParticle, prevParticle, second_prevParticle;
	double dt, dt_estimate, mag, curTime;
	VECTOR3 vel;
	double cell_volume;
	int save_counter = 0; // only if the save_counter == 0 will store the point

	// the first particle
	istat = m_pField->at_phys(seedInfo.fromCell, seedInfo.phyCoord, seedInfo, m_fCurrentTime, vel);
	if(istat == OUT_OF_BOUND)  {
		return OUT_OF_BOUND;
	}
	// Do not check critical point for MPAS Ocean data.
	// if((fabs(vel[0]) < m_fStationaryCutoff) && (fabs(vel[1]) < m_fStationaryCutoff) && (fabs(vel[2]) < m_fStationaryCutoff)) {
	// 	return CRITICAL_POINT;
	// }

	thisParticle = seedInfo;

	seedTrace.push_back(new VECTOR3(seedInfo.phyCoord));
	curTime = m_fCurrentTime;
	count++;

	// get the initial stepsize
	// this is a bug I think ...
	if(m_pField->GetCellType() == VORONOI) {
		dt = dt_estimate = m_fInitStepSize;
	}
	else {
		cell_volume = m_pField->volume_of_cell(seedInfo.inCell);
		mag = vel.GetMag();
// 		printf(" **** volume = %f mag = %f  ....", cell_volume, mag);
		if(fabs(mag) < 1.0e-6f)
			dt_estimate = 1.0e-5f;
		else
			dt_estimate = pow(cell_volume, (double)0.3333333f) / mag;
		dt = m_fInitStepSize * dt_estimate;
	}
			//	printf(" **** dt = %f  ....", dt);

#ifdef DEBUG
	fprintf(fDebugOut, "****************new particle*****************\n");
	fprintf(fDebugOut, "seed: %f, %f, %f with step size %f\n", seedInfo.phyCoord[0], seedInfo.phyCoord[1], seedInfo.phyCoord[2], dt);
#endif

	// start to advect
	while(count < m_nMaxsize)
	{
		second_prevParticle = prevParticle;
		prevParticle = thisParticle;

		if(integ_ord == SECOND)
			istat = runge_kutta2(time_dir, time_dep, thisParticle, &curTime, dt);
		else if(integ_ord == MPASO_EULER) {
			dt = m_fInitStepSize;
			istat = MPASO_euler(time_dir, time_dep, thisParticle, &curTime, dt);
		}
		else if(integ_ord == MPASO_FOURTH) {
			dt = m_fInitStepSize;
			istat = MPASO_rk4(time_dir, time_dep, thisParticle, &curTime, dt);
		}
		else
			istat = runge_kutta4(time_dir, time_dep, thisParticle, &curTime, dt);

		if(istat == OUT_OF_BOUND) {			// out of boundary
		    seedTrace.push_back(new VECTOR3(thisParticle.phyCoord));
			return OUT_OF_BOUND;
		}
		m_pField->resolve_cell(thisParticle.fromCell, thisParticle.phyCoord, thisParticle);
		// Only store points every 1000 step

		// if(save_counter == 0) {
			seedTrace.push_back(new VECTOR3(thisParticle.phyCoord));
			count++;	
		// }
		// save_counter++;
		// if(save_counter == 1000) save_counter = 0;

		// MPASO_EULER resets dt to a fixed value at the top of the loop,
		// so adapt_step is irrelevant for that integrator and is skipped.
		if(count > 2 && integ_ord != MPASO_EULER && integ_ord != MPASO_FOURTH)
		{
			adapt_step(second_prevParticle.phyCoord, prevParticle.phyCoord, thisParticle.phyCoord, dt_estimate, &dt);
		}
	}

	return OKAY;
}


int vtCStreamLine::computeFieldLine( TIME_DIR time_dir,
									 INTEG_ORD integ_ord,
									 TIME_DEP time_dep,
									 vtListSeedTrace& seedTrace,
									 PointInfo& seedInfo,
									 int fromCell,
									 vector<int>& toCells)
{
	int count = 0, istat;
	PointInfo thisParticle, prevParticle, second_prevParticle;
	double dt, dt_estimate, mag, curTime;
	VECTOR3 vel;
	double cell_volume;
	int stepCount = 0;       // actual integration steps taken
	bool lastStored = true;  // seed point counts as "stored" for OUT_OF_BOUND logic

	// the first particle
	if(fromCell >= 0)
		seedInfo.fromCell = fromCell;
	istat = m_pField->at_phys(seedInfo.fromCell, seedInfo.phyCoord, seedInfo, m_fCurrentTime, vel);
	if(istat == OUT_OF_BOUND)  {
		toCells.push_back(-1);  // keep toCells in sync with sl_list
		m_lastFieldLineStepCount = 0;
		return OUT_OF_BOUND;
	}

	thisParticle = seedInfo;

	seedTrace.push_back(new VECTOR3(seedInfo.phyCoord));
	curTime = m_fCurrentTime;
	count++;

	// get the initial stepsize
	// this is a bug I think ...
	if(m_pField->GetCellType() == VORONOI) {
		dt = dt_estimate = m_fInitStepSize;
	}
	else {
		cell_volume = m_pField->volume_of_cell(seedInfo.inCell);
		mag = vel.GetMag();
		if(fabs(mag) < 1.0e-6f)
			dt_estimate = 1.0e-5f;
		else
			dt_estimate = pow(cell_volume, (double)0.3333333f) / mag;
		dt = m_fInitStepSize * dt_estimate;
	}

	// start to advect
	while(count < m_nMaxsize)
	{
		second_prevParticle = prevParticle;
		prevParticle = thisParticle;

		if(integ_ord == SECOND)
			istat = runge_kutta2(time_dir, time_dep, thisParticle, &curTime, dt);
		else if(integ_ord == MPASO_EULER) {
			dt = m_fInitStepSize;
			istat = MPASO_euler(time_dir, time_dep, thisParticle, &curTime, dt);
		}
		else if(integ_ord == MPASO_FOURTH) {
			dt = m_fInitStepSize;
			istat = MPASO_rk4(time_dir, time_dep, thisParticle, &curTime, dt);
		}
		else
			istat = runge_kutta4(time_dir, time_dep, thisParticle, &curTime, dt);

		if(istat == OUT_OF_BOUND) {
			toCells.push_back(thisParticle.fromCell);
			// Ensure the last valid position (prevParticle) is in the trace for re-injection
			if (!lastStored)
				seedTrace.push_back(new VECTOR3(prevParticle.phyCoord));
			m_lastFieldLineStepCount = stepCount;
			return OUT_OF_BOUND;
		}
		m_pField->resolve_cell(thisParticle.fromCell, thisParticle.phyCoord, thisParticle);

		stepCount++;
		lastStored = (stepCount % m_nSaveInterval == 0);
		if (lastStored)
			seedTrace.push_back(new VECTOR3(thisParticle.phyCoord));
		count++;

		if(count > 2)
		{
			double tempt;
			tempt = dt;
			adapt_step(second_prevParticle.phyCoord, prevParticle.phyCoord, thisParticle.phyCoord, dt_estimate, &dt);
		}
	}

	// Always store the final position if not already stored (needed for re-injection)
	if (!lastStored)
		seedTrace.push_back(new VECTOR3(thisParticle.phyCoord));
	m_lastFieldLineStepCount = stepCount;
	toCells.push_back(thisParticle.fromCell);

	return OKAY;
}

// streamline advects as far as possible till the boundary or terminates at critical points
// only two cases will happen: OUT_OF_BOUNDARY & CRITICAL_POINT
int vtCStreamLine::executeInfiniteAdvection(TIME_DIR time_dir,
											TIME_DEP time_dep,
											vtListSeedTrace& seedTrace,
											double& totalStepsize,
											vector<double>* vStepsize)
{
	int istat;
	vtParticleInfo* thisSeed;
	PointInfo thisParticle, prevParticle, second_prevParticle, seedInfo;
	vtListParticleIter sIter;
	double dt, dt_estimate, cell_volume, mag, curTime;
	VECTOR3 vel;
	INTEG_ORD integ_ord;

	// initialize
	integ_ord = m_integrationOrder;
	sIter = m_lSeeds.begin();
	thisSeed = *sIter;
	seedInfo = thisSeed->m_pointInfo;
	seedTrace.clear();
	totalStepsize = 0.0;

	// the first particle
	seedTrace.push_back(new VECTOR3(seedInfo.phyCoord));
	istat = m_pField->at_phys(seedInfo.fromCell, seedInfo.phyCoord, seedInfo, m_fCurrentTime, vel);
	if((fabs(vel[0]) < m_fStationaryCutoff) && (fabs(vel[1]) < m_fStationaryCutoff) && (fabs(vel[2]) < m_fStationaryCutoff))
		return CRITICAL_POINT;

	thisParticle = seedInfo;
	curTime = m_fCurrentTime;

	// get the initial stepsize
	if(m_pField->GetCellType() == VORONOI) {
		dt = dt_estimate = m_fInitStepSize;
	}
	else {
		cell_volume = m_pField->volume_of_cell(seedInfo.inCell);
		mag = vel.GetMag();
		dt_estimate = pow(cell_volume, (double)0.3333333f) / mag;
		dt = m_fInitStepSize * dt_estimate;
	}

#ifdef DEBUG
	fprintf(fDebugOut, "****************new particle*****************\n");
	fprintf(fDebugOut, "seed: %f, %f, %f with step size %f\n", seedInfo.phyCoord[0], seedInfo.phyCoord[1], seedInfo.phyCoord[2], dt);
#endif

	// start to advect
	while(true)
	{
		second_prevParticle = prevParticle;
		prevParticle = thisParticle;

		if(integ_ord == SECOND)
			istat = runge_kutta2(time_dir, time_dep, thisParticle, &curTime, dt);
		else
			istat = runge_kutta4(time_dir, time_dep, thisParticle, &curTime, dt);

		if(istat == OUT_OF_BOUND)			// out of boundary
		{
			// find the boundary intersection point
			VECTOR3 intersectP, startP, endP;
			double oldStepsize, stepSize; 
			oldStepsize = stepSize = dt;
			startP = prevParticle.phyCoord;
			endP = thisParticle.phyCoord;
			m_pField->BoundaryIntersection(intersectP, startP, endP, &stepSize, oldStepsize);
			totalStepsize += stepSize;
			seedTrace.push_back(new VECTOR3(intersectP));
			if(vStepsize != NULL)
				vStepsize->push_back(stepSize);
			return OUT_OF_BOUND;
		}

		// find the loop
		list<VECTOR3*>::iterator searchIter;
		searchIter = seedTrace.end();
		searchIter--;
		for(; searchIter != seedTrace.begin(); searchIter--)
		{
			if(thisParticle.phyCoord == **searchIter)			// loop
				return CRITICAL_POINT;
		}

		m_pField->at_phys(thisParticle.fromCell, thisParticle.phyCoord, thisParticle, m_fCurrentTime, vel);
		seedTrace.push_back(new VECTOR3(thisParticle.phyCoord));
		if(vStepsize != NULL)
			vStepsize->push_back(dt);
		totalStepsize += dt;

#ifdef DEBUG
		fprintf(fDebugOut, "****************advected particle*****************\n");
		fprintf(fDebugOut, "pos = (%f, %f, %f), vel = (%f, %f, %f) with step size %f, total step size %f\n", thisParticle.phyCoord[0], thisParticle.phyCoord[1], thisParticle.phyCoord[2], vel[0], vel[1], vel[2], dt, totalStepsize);
#endif

		if(totalStepsize > 4000.0)
		{
			printf("curl\n");
			vel.Set(0.0, 0.0, 0.0);
		}

		if((fabs(vel[0]) < m_fStationaryCutoff) && (fabs(vel[1]) < m_fStationaryCutoff) && (fabs(vel[2]) < m_fStationaryCutoff))
			return CRITICAL_POINT;
		if((int)seedTrace.size() > 2)
			adapt_step(second_prevParticle.phyCoord, prevParticle.phyCoord, thisParticle.phyCoord, dt_estimate, &dt);
	}

	return OUT_OF_BOUND;
}

int vtCStreamLine::AdvectOneStep(TIME_DIR time_dir,
								 INTEG_ORD integ_ord,
								 TIME_DEP time_dep,
								 PointInfo& seedInfo,
								 VECTOR3& finalP)
{
	int istat;
	PointInfo thisParticle;
	VECTOR3 vel;
	double curTime = m_fCurrentTime;
	double dt = m_fInitStepSize;

	finalP.Set(-1, -1, -1);
	thisParticle = seedInfo;
	
	// the first particle
	istat = m_pField->at_phys(seedInfo.fromCell, seedInfo.phyCoord, seedInfo, m_fCurrentTime, vel);
	if(istat == OUT_OF_BOUND)
		return OUT_OF_BOUND;
	if((fabs(vel[0]) < m_fStationaryCutoff) && (fabs(vel[1]) < m_fStationaryCutoff) && (fabs(vel[2]) < m_fStationaryCutoff))
		return CRITICAL_POINT;
	if(integ_ord == SECOND)
		istat = runge_kutta2(time_dir, time_dep, thisParticle, &curTime, dt);
	else
		istat = runge_kutta4(time_dir, time_dep, thisParticle, &curTime, dt);

	if(istat == OUT_OF_BOUND)			// out of boundary
			return OUT_OF_BOUND;
	m_pField->at_phys(thisParticle.fromCell, thisParticle.phyCoord, thisParticle, m_fCurrentTime, vel);
	if((fabs(vel[0]) < m_fStationaryCutoff) && (fabs(vel[1]) < m_fStationaryCutoff) && (fabs(vel[2]) < m_fStationaryCutoff))
		return CRITICAL_POINT;
	else
		finalP = thisParticle.phyCoord;
	return OKAY;
}
