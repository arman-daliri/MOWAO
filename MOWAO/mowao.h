#ifndef __MOWAO__
#define __MOWAO__

typedef void (*problem_func)(double *x, int x_len, double *f);

struct Particle {
	double *h; /* decision variables */
	double *velocity;
	double *pl; /* objectives */
	int evaporated;
	int numh; /* number of particles */
};

struct Repository {
	struct Particle *repo;
	double *dist;
	double *max;
	int end;
};

struct Population {
	struct Particle *pop;
	int end;
};

struct MOWAO {
	int nobj; /* objectives */
	int nvar;
	double *lb, *ub;
	int nrepo; /* repo size */
	int hpop;
	int maxiter;

	/* parameters */
	double bond_radius;
	double punch;
	double evaporate;
	double coef;
	double *vlb, *vub;

	problem_func f;
	//void (*f)(double *x, int x_len, double *pl);
	struct Population pop;
	struct Repository rep;
};

double randlim(double min, double max);
int dominates(double *p1, double *p2, int nobj);
void check_boundary(double *h, int n, double *ub, double *lb);
double distance(double *p1, double *p2, int len, double *max);

struct MOWAO *mowao_new(void);
void mowao_alloc(struct MOWAO *mw);
void mowao_init(struct MOWAO *mw);
void mowao_run(struct MOWAO *mw, int log);
void mowao_clean(struct MOWAO *mw);
void mowao_free(struct MOWAO *mw);

void particle_alloc(struct MOWAO *mw, struct Particle *p);
void particle_init(struct MOWAO *mw, struct Particle *p);
void particle_clean(struct Particle *p);
void particle_copy(struct MOWAO *mw, struct Particle *dst,
		   struct Particle *src);
void particle_punch(struct MOWAO *mw, struct Particle *p);
void particle_evaporate(struct MOWAO *mw, struct Particle *p);

void population_init(struct MOWAO *mw);
void population_clean(struct MOWAO *mw);
void population_update(struct MOWAO *mw);
void population_fill(struct MOWAO *mw);

void repository_init(struct MOWAO *mw);
void repository_clean(struct MOWAO *mw);
void repository_check(struct MOWAO *mw);
void repository_add(struct MOWAO *mw, struct Particle *p);
void repository_update(struct MOWAO *mw, struct Particle *p);
double repository_dist(struct MOWAO *mw, struct Particle *p);

#endif /* __MOWAO__ */
