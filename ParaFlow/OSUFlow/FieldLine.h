/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 Vector Field Data
//
///////////////////////////////////////////////////////////////////////////////


#ifndef _VECTOR_FIELD_LINE_H_
#define _VECTOR_FIELD_LINE_H_

#include "Field.h"

//Silence an annoying and unnecessary compiler warning
#pragma warning(disable : 4251 4100 4244)

//////////////////////////////////////////////////////////////////////////
// definition
//////////////////////////////////////////////////////////////////////////
#define MAX_LENGTH 1000000
const double STREAM_ACCURACY = EPS;
#define INFINITE_LINE -1				// advect fieldlines as far as possible

enum INTEG_ORD{ SECOND = 2, FOURTH = 4, MPASO_EULER = 11, MPASO_FOURTH = 14};		// integration order
enum TIME_DIR{ BACKWARD = -1, FORWARD = 1};		// advection direction
enum TIME_DEP{ STEADY=0,UNSTEADY=1 };	
enum TRACE_DIR{OFF=0, BACKWARD_DIR=1, FORWARD_DIR=2, BACKWARD_AND_FORWARD=3};
enum ADVECT_STATUS{NONE = -3, OUT_OF_TIME = -2, OUT_OF_BOUND = -1, CRITICAL_POINT = 0, OKAY = 1};

//////////////////////////////////////////////////////////////////////////
// information about particles
//////////////////////////////////////////////////////////////////////////
class vtParticleInfo
{
public:
	PointInfo m_pointInfo;		// basic information about this particle
	double m_fStartTime;			// start time
	double m_fCurrentTime;		// current time inside the active block
	int m_cachedLowT;			// per-particle local time-interval hint; -1 = invalid
	int itsValidFlag;			// whether this particle is valid or not
	int itsNumStepsAlive;		// number of steps alive
	int ptId;					// particle ID

public:
	vtParticleInfo(void)
	{
		m_fStartTime = 0.0;
		m_fCurrentTime = 0.0;
		m_cachedLowT = -1;
		itsValidFlag = 0;
		itsNumStepsAlive = 0;
		ptId = -1;
	}

	vtParticleInfo(vtParticleInfo* source)
	{
		m_pointInfo.Set(source->m_pointInfo);
		m_fStartTime = source->m_fStartTime;
		m_fCurrentTime = source->m_fCurrentTime;
		m_cachedLowT = source->m_cachedLowT;
		itsValidFlag = source->itsValidFlag;
		itsNumStepsAlive = 0;
		ptId = -1;
	}

	vtParticleInfo(vtParticleInfo& source)
	{
		m_pointInfo.Set(source.m_pointInfo);
		m_fStartTime = source.m_fStartTime;
		m_fCurrentTime = source.m_fCurrentTime;
		m_cachedLowT = source.m_cachedLowT;
		itsValidFlag = source.itsValidFlag;
		ptId = source.ptId;
		itsNumStepsAlive = 0;
	}
};

typedef list<vtParticleInfo*> vtListParticle;
typedef list<vtParticleInfo*>::iterator vtListParticleIter;
typedef list<VECTOR3*> vtListSeedTrace;      // positions
typedef list<VECTOR4*> vtListTimeSeedTrace;  // positions and times 

//////////////////////////////////////////////////////////////////////////
// base class for all field lines, and the structure is like:
//                                    vtCFieldLine
//	                                  /  \
//			              vtCStreamLine   vtCTimeVaryingFieldLine
//               	          	|              /   \
//                      vtCPathline  vtCTimeLine   vtCStreakLine
//////////////////////////////////////////////////////////////////////////
class vtCFieldLine
{
protected:
	int m_nNumSeeds;		// number of seeds
	INTEG_ORD m_integrationOrder;	// integration order
	TIME_DIR m_timeDir;	        // advection direction
	TIME_DEP m_itsTimeDep;
	double m_fInitTime;
	double m_fStepTime;
	double m_fInitStepSize;	 // initial advection step size of particle
	double m_fDurationTime;
	double m_fLowerAngleAccuracy;	// for adaptive stepsize 
	double m_fUpperAngleAccuracy;
	double m_fMaxStepSize;	        // maximal advection stepsize
	int m_nMaxsize;		// maximal number of particles this line advects
	int m_nSaveInterval;    // store every N-th integration step (1 = store all)
	int m_lastFieldLineStepCount; // actual steps taken in the last computeFieldLine/advectParticle call
	std::vector<int> m_traceActualStepCounts; // per-trace step counts accumulated by computeStreamLine/computePathLine
	vtListParticle m_lSeeds;	// list of seeds
	list<int64_t> m_lSeedIds;	// list of seed ids
	CVectorField* m_pField;	        // vector field
	double m_fStationaryCutoff;	// cutoff value for critical points

public:
	vtCFieldLine(CVectorField* pField);
	virtual ~vtCFieldLine(void);

	void setSeedPoints(VECTOR3* points, int numPoints, double t, int64_t *seedIds = NULL); 
	void setSeedPoints(VECTOR3* points, int numPoints, double* tarray);
	void setSeedPoints(VECTOR4* points, int numPoints, int64_t *seedIds = NULL);
	void setMaxPoints(int val) { m_nMaxsize = val; }
	void setSaveInterval(int n) { m_nSaveInterval = (n < 1) ? 1 : n; }
	const std::vector<int>& getTraceActualStepCounts() const { return m_traceActualStepCounts; }
	void setIntegrationOrder(INTEG_ORD ord) { m_integrationOrder = ord; }
	int  getMaxPoints(void){ return m_nMaxsize; }
	INTEG_ORD getIntegrationOrder(void){ return m_integrationOrder; }
	void SetMaxStepSize(double stepsize) {m_fMaxStepSize = stepsize;}
	double GetMaxStepSize(void) {return m_fMaxStepSize;}
	void SetInitStepSize(double initStep) { m_fInitStepSize = initStep; }
	double GetInitStepSize(void) { return m_fInitStepSize; }
	void SetLowerUpperAngle(double lowerAngle, double upperAngle) {m_fLowerAngleAccuracy = lowerAngle; m_fUpperAngleAccuracy = upperAngle;}
	void SetStationaryCutoff(double cutoff) {m_fStationaryCutoff = cutoff;}

protected:
	void releaseSeedMemory(void);
	int euler_cauchy(TIME_DIR, TIME_DEP,double*, double);
	int runge_kutta4(TIME_DIR, TIME_DEP, PointInfo&, double*, double, int* cachedLowT = NULL);
	int runge_kutta2(TIME_DIR, TIME_DEP, PointInfo&, double*, double, int* cachedLowT = NULL);
	int  MPASO_euler(TIME_DIR, TIME_DEP, PointInfo&, double*, double, int* cachedLowT = NULL);
	int  MPASO_rk4(TIME_DIR, TIME_DEP, PointInfo&, double*, double, int* cachedLowT = NULL);
	bool geodesic_step(TIME_DIR, VECTOR3, double, VECTOR3, double, double, VECTOR3&, double&);
	int adapt_step(const VECTOR3& p2, const VECTOR3& p1, const VECTOR3& p0, double dt_estimate,double* dt);
};

//////////////////////////////////////////////////////////////////////////
// class declaratin of timevaryingfieldline
//////////////////////////////////////////////////////////////////////////
class vtCTimeVaryingFieldLine : public vtCFieldLine
{
protected:
	int m_itsTimeAdaptionFlag;
	int m_itsMaxParticleLife;   // how long the particles be alive
	int m_itsMapWithTimeFlag;
	int m_itsNeedResetFlag;
	int m_itsWrapTimeFlag;
	double m_itsTimeInc;
	list<vtParticleInfo*> m_itsParticles;
public:
	vtCTimeVaryingFieldLine(CVectorField* pField);
	~vtCTimeVaryingFieldLine(void);

	void SetTimeDir(TIME_DIR dir) { m_timeDir = dir; }

	void SetInjectionTime(double init, double step, double duration, TIME_DIR dir = FORWARD)
	{
		m_fInitTime = init;
		m_fStepTime = step;
		m_fDurationTime = duration;
		m_timeDir = dir;
	}

	void SetTimeAdaptionMode(int onoff)
	{
		m_itsTimeAdaptionFlag = onoff;
	}

	virtual void setParticleLife(int steps);
	virtual int getParticleLife(void);
	virtual void setTimeMapping(int enabled);
	virtual int getTimeMapping(void);
	virtual void killAllParticles(void) { m_itsNeedResetFlag = 1; }

protected:
	// code shared by all time varying field lines here
	void releaseParticleMemory(void);
	int advectParticle( INTEG_ORD int_order, 
			    vtParticleInfo& initialPoint,
			    double initialTime,
			    vtParticleInfo& finalPoint,
			    double finalTime);
	int advectParticle( INTEG_ORD int_order, 
			    vtParticleInfo& initialPoint,
			    double initialTime,
			    double finalTime,
			    vtListSeedTrace& seedTrace);
	int advectParticle( INTEG_ORD int_order, 
			    vtParticleInfo& initialPoint,
			    double initialTime,
			    double finalTime,
			    vtListTimeSeedTrace& seedTrace,
			    vector<int>* toCells = nullptr);
};

//////////////////////////////////////////////////////////////////////////
// class declaration of pathline
//////////////////////////////////////////////////////////////////////////
typedef struct vtPathlineParticle
{
	VECTOR3 pos;
	int ptId;
        double time; 
}vtPathlineParticle;

class vtCPathLine : public vtCTimeVaryingFieldLine
{
public:
	vtCPathLine(CVectorField* pField);
	~vtCPathLine(void);

	void execute(list<vtListTimeSeedTrace*>& listSeedTraces, list<int64_t> *listSeedIds = NULL);
	void execute(list<vtListTimeSeedTrace*>& listSeedTraces, vector<int>& fromCells, vector<int>& toCells);
	void execute(list<vtPathlineParticle*>& listSeedTraces);

protected:
	// code specific to pathline
	void computePathLine(list<vtListTimeSeedTrace*>&, list<int64_t> *listSeedIds = NULL);
	void computePathLine(list<vtListTimeSeedTrace*>&, vector<int>& fromCells, vector<int>& toCells);
	void computePathLine(list<vtPathlineParticle*>&);
};

//////////////////////////////////////////////////////////////////////////
// class declaration of streakline
// inject particles from a point
//////////////////////////////////////////////////////////////////////////
typedef struct vtStreakParticle
{
	PointInfo itsPoint;			// basic information
	double itsTime;				// start time
	int traceId;				// which tracing path this particle belongs to, used to track advection
}vtStreakParticle;

typedef vector<vtStreakParticle*> vtListStreakParticle;	    // one pathline trace from a particle
typedef vector<vtListStreakParticle*> vtStreakTraces;	    // all pathline traces from all particles released
typedef vector<vtStreakParticle*>::iterator vtStreakParticleIter;
typedef vector<vtListStreakParticle*>::iterator vtStreakTracesIter; 

class vtCStreakLine : public vtCTimeVaryingFieldLine
{
public:
	vtCStreakLine(CVectorField* pField);
	~vtCStreakLine(void);
	void execute(const void* userData, vtStreakTraces& listSeedTraces);

protected:
	// code specific to streakline
	void computeStreakLine(const void*, vtStreakTraces&);
	void advectOldParticles(vtListParticleIter start, 
				vtListParticleIter end, 
				vtStreakTraces& listSeedTraces,
				double initialTime, double finalTime,
				vector<vtListParticleIter>& deadList);
private:
	int nHowManyTraces;
};

//////////////////////////////////////////////////////////////////////////
// class declaration of timeline
// release a set of particles, injected periodically
//////////////////////////////////////////////////////////////////////////
class vtCTimeLine : public vtCTimeVaryingFieldLine
{
public:
	vtCTimeLine(CVectorField* pField);
	~vtCTimeLine(void);

	void execute(const void* userData, vtListStreakParticle& listSeedTraces);
	void setTimeDelay(int delay);
	int getTimeDelay(void);

protected:
	// code specific to timeline
	void computeTimeLine(const void*, vtListStreakParticle&);
	void advectOldParticles(vtListParticleIter start, 
				vtListParticleIter end, 
				vtListStreakParticle& listSeedTraces,
				double initialTime, double finalTime,
				vector<vtListParticleIter>& deadList);
	int m_itsTimeDelay;
	int numTillRelease;
};

//////////////////////////////////////////////////////////////////////////
// class declaration of streamline
//////////////////////////////////////////////////////////////////////////
class vtCStreamLine : public vtCFieldLine
{
public:
	vtCStreamLine(CVectorField* pField);
	~vtCStreamLine(void);

	void execute(const void* userData, list<vtListSeedTrace*>& listSeedTraces,
				list<int64_t> *listSeedIds = NULL);
	void execute(const void* userData, list<vtListSeedTrace*>& listSeedTraces,
				vector<int>& fromCells, vector<int>& toCells, MPASOReader* mpaso_reader, list<int64_t> *listSeedIds = NULL);
	void setForwardTracing(int enabled);
	void setBackwardTracing(int enabled);
	int  getForwardTracing(void);
	int  getBackwardTracing(void);
	int executeInfiniteAdvection(TIME_DIR, TIME_DEP, vtListSeedTrace&, double&, vector<double>*);
	int AdvectOneStep(TIME_DIR, INTEG_ORD, TIME_DEP, PointInfo&, VECTOR3&);

protected:
	void computeStreamLine(const void* userData, list<vtListSeedTrace*>& listSeedTraces, list<int64_t> *listSeedIds = NULL);
	void computeStreamLine(const void* userData, list<vtListSeedTrace*>& listSeedTraces, vector<int>& fromCells, vector<int>& toCells, MPASOReader* mpaso_reader, list<int64_t> *listSeedIds = NULL);
	int computeFieldLine(TIME_DIR, INTEG_ORD, TIME_DEP, vtListSeedTrace&, PointInfo&);
	int computeFieldLine(TIME_DIR, INTEG_ORD, TIME_DEP, vtListSeedTrace&, PointInfo&, int, vector<int>&);

	TRACE_DIR m_itsTraceDir;
	double m_fPsuedoTime;
	double m_fCurrentTime;
};

#endif
