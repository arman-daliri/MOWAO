#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "mowao.h"

double
randlim(double min, double max)
{
	double range = (max - min);
	double div = RAND_MAX / range;
	return min + (rand() / div);
}

int
dominates(double *p1, double *p2, int nobj)
{
	int i, s = 0;
	for (i = 0; i < nobj; i++) {
		if (p1[i] > p2[i])
			return 0;
		else if (p1[i] != p2[i])
			s = 1;
	}
	return s;
}

void
check_boundary(double *h, int len, double *ub, double *lb)
{
	int i;

	for (i = 0; i < len; i++) {
		if (h[i] > ub[i])
			h[i] = ub[i];
		if (h[i] < lb[i])
			h[i] = lb[i];
	}
}

double
distance(double *p1, double *p2, int len, double *max)
{
	int i;
	double d = 0;

	for (i = 0; i < len; i++)
		d += pow((p1[i] - p2[i]) / (max[i] + DBL_MIN), 2);
	return sqrt(d);
}

struct MOWAO *
mowao_new(void)
{
	struct MOWAO *mw = calloc(1, sizeof(struct MOWAO));
	return mw;
}

void
mowao_alloc(struct MOWAO *mw)
{
	mw->lb = calloc(mw->ndec, sizeof(double));
	mw->ub = calloc(mw->ndec, sizeof(double));
	mw->vlb = calloc(mw->ndec, sizeof(double));
	mw->vub = calloc(mw->ndec, sizeof(double));
}

void
mowao_init(struct MOWAO *mw)
{
	srand(time(NULL));
	population_init(mw);
	repository_init(mw);
	repository_update(mw, mw->pop.pop);
}

void
mowao_run(struct MOWAO *mw, int log)
{
	int i;

	for (i = 0; i < mw->maxiter; i++) {
		if (log) {
			printf("Iteration: %d, Population: %d, Repository: %d\n",
			       i + 1, mw->pop.end, mw->rep.end);
		}
		population_update(mw);
		repository_update(mw, mw->pop.pop);
	}
}

void
mowao_clean(struct MOWAO *mw)
{
	population_clean(mw);
	repository_clean(mw);
	mowao_clean_memebers(mw);
	mowao_clean_bounds(mw);
}

void
mowao_clean_memebers(struct MOWAO *mw)
{
	int i;

	for (i = 0; i > mw->ndec; i++) {
		mw->vlb[i] = 0;
		mw->vub[i] = 0;
		mw->lb[i] = 0;
		mw->ub[i] = 0;
	}
	mw->nobj = 0;
	mw->ndec = 0;
	mw->nrepo = 0;
	mw->npop = 0;
	mw->maxiter = 0;
	mw->bond_radius = 0;
	mw->push = 0;
	mw->evaporate = 0;
	mw->coef = 0;
}

void
mowao_clean_bounds(struct MOWAO *mw)
{
	free(mw->lb);
	mw->lb = NULL;
	free(mw->ub);
	mw->ub = NULL;
	free(mw->vlb);
	mw->vlb = NULL;
	free(mw->vub);
	mw->vub = NULL;
}

void
mowao_free(struct MOWAO *mw)
{
	free(mw);
}

void
particle_alloc(struct MOWAO *mw, struct Particle *p)
{
	p->dec = calloc(mw->ndec, sizeof(double));
	p->velocity = calloc(mw->ndec, sizeof(double));
	p->f = calloc(mw->nobj, sizeof(double));
}

void
particle_init(struct MOWAO *mw, struct Particle *p)
{
	int i;

	for (i = 0; i < mw->ndec; i++) {
		p->dec[i] = randlim(mw->lb[i], mw->ub[i]);
		p->velocity[i] = randlim(mw->vlb[i], mw->vub[i]);
	}
	mw->f(p->dec, mw->ndec, p->f, mw->cb_arg);
	p->num_bond = 1;
	p->evaporated = 0;
}

void
particle_clean(struct Particle *p)
{
	free(p->dec);
	p->dec = NULL;
	free(p->velocity);
	p->velocity = NULL;
	free(p->f);
	p->f = NULL;
}

void
particle_copy(struct MOWAO *mw, struct Particle *dst, struct Particle *src)
{
	memcpy(dst->dec, src->dec, mw->ndec * sizeof(double));
	memcpy(dst->velocity, src->velocity, mw->ndec * sizeof(double));
	memcpy(dst->f, src->f, mw->nobj * sizeof(double));

	dst->num_bond = src->num_bond;
}

void
particle_push(struct MOWAO *mw, struct Particle *p)
{
	int i, ind = 0;
	double rnd1, rnd2;

	ind = randlim(0, mw->rep.end);
	rnd1 = randlim(0, 1);
	rnd2 = randlim(0, 1);
	for (i = 0; i < mw->ndec; i++) {
		p->velocity[i] = rnd1 * p->velocity[i];
		p->velocity[i] += rnd2 * (mw->rep.repo[ind].dec[i] - p->dec[i]);
	}
}

void
particle_evaporate(struct MOWAO *mw, struct Particle *p)
{
	int i;

	for (; p->num_bond > 1; p->num_bond--) {
		mw->pop.pop[mw->pop.end].num_bond = 1;
		mw->pop.pop[mw->pop.end].evaporated = 1;
		for (i = 0; i < mw->ndec; i++) {
			mw->pop.pop[mw->pop.end].dec[i] = p->dec[i];
			mw->pop.pop[mw->pop.end].velocity[i] =
				randlim(mw->vlb[i], mw->vub[i]);
			mw->pop.pop[mw->pop.end].velocity[i] *= mw->coef;
		}
		mw->pop.end++;
	}
	p->evaporated = 1;
	for (i = 0; i < mw->ndec; i++) {
		p->velocity[i] = randlim(mw->vlb[i], mw->vub[i]);
		p->velocity[i] *= mw->coef;
	}
}

void
population_init(struct MOWAO *mw)
{
	int i;
	mw->pop.pop = calloc(mw->npop, sizeof(struct Particle));
	mw->pop.end = 0;
	for (i = 0; i < mw->npop; i++)
		particle_alloc(mw, &(mw->pop.pop[i]));
	population_fill(mw);
}

void
population_clean(struct MOWAO *mw)
{
	int i;

	for (i = 0; i < mw->npop; i++)
		particle_clean(&(mw->pop.pop[i]));
	free(mw->pop.pop);
	mw->pop.pop = NULL;
}

void
population_fill(struct MOWAO *mw)
{
	for (; mw->pop.end < mw->npop; mw->pop.end++)
		particle_init(mw, &(mw->pop.pop[mw->pop.end]));
}

void
population_update(struct MOWAO *mw)
{
	int i, j, k;
	struct Particle *p = mw->pop.pop;

	for (i = 0; i < mw->pop.end; i++) {
		for (j = i + 1; j < mw->pop.end; j++) {
			if (distance(p[i].f, p[j].f, mw->nobj, mw->rep.max) <=
			    mw->bond_radius) {
				for (k = 0; k < mw->ndec; k++)
					p[i].velocity[k] += p[j].velocity[k];

				p[i].num_bond++;
				particle_copy(mw, &p[j], &p[mw->pop.end - 1]);
				mw->pop.end--;
				j--;
			}
		}
	}

	double perc = (double)mw->pop.end / mw->npop * 100;
	if (perc < 70) {
		mw->evaporate += 0.1;
	} else if (perc > 90) {
		mw->evaporate -= 0.1;
	}
	if (mw->evaporate > 1) {
		mw->evaporate = 1;
	} else if (mw->evaporate < 0) {
		mw->evaporate = 0;
	}

	for (i = 0; i < mw->pop.end; i++) {
		if (p[i].num_bond > 1 && randlim(0, 1) <= mw->evaporate)
			particle_evaporate(mw, &p[i]);
		else if (randlim(0, 1) <= mw->push)
			particle_push(mw, &p[i]);
		if (p[i].evaporated == 0)
			check_boundary(p[i].velocity, mw->ndec, mw->vub,
				       mw->vlb);
		else
			p[i].evaporated = 0;
		for (j = 0; j < mw->ndec; j++) {
			p[i].dec[j] += p[i].velocity[j];
		}
		check_boundary(p[i].dec, mw->ndec, mw->ub, mw->lb);
		mw->f(p[i].dec, mw->ndec, p[i].f, mw->cb_arg);
	}
}

void
repository_init(struct MOWAO *mw)
{
	int i;

	mw->rep.repo = calloc(mw->nrepo, sizeof(struct Particle));
	mw->rep.dist = calloc(mw->nrepo, sizeof(double));
	mw->rep.max = calloc(mw->nobj, sizeof(double));

	for (i = 0; i < mw->nrepo; i++)
		particle_alloc(mw, &(mw->rep.repo[i]));
	mw->rep.end = 0;
}

void
repository_clean(struct MOWAO *mw)
{
	int i;

	for (i = 0; i < mw->nrepo; i++)
		particle_clean(&(mw->rep.repo[i]));
	free(mw->rep.repo);
	mw->rep.repo = NULL;
	free(mw->rep.dist);
	mw->rep.dist = NULL;
	free(mw->rep.max);
	mw->rep.max = NULL;
}

void
repository_check(struct MOWAO *mw)
{
	int i, j;

	for (i = 0; i < mw->rep.end; i++) {
		for (j = 0; j < mw->rep.end; j++) {
			if (i == j)
				continue;
			if (dominates(mw->rep.repo[i].f, mw->rep.repo[j].f,
				      mw->nobj)) {
				particle_copy(mw, &(mw->rep.repo[j]),
					      &(mw->rep.repo[mw->rep.end - 1]));
				mw->rep.end--;
				j--;
				continue;
			}
			if (distance(mw->rep.repo[i].f, mw->rep.repo[j].f,
				     mw->nobj,
				     mw->rep.max) <= mw->bond_radius) {
				mw->rep.repo[i].num_bond++;
				particle_copy(mw, &(mw->rep.repo[j]),
					      &(mw->rep.repo[mw->rep.end - 1]));
				mw->rep.end--;
				j--;
			}
		}
	}
}

void
repository_add(struct MOWAO *mw, struct Particle *p)
{
	int i, ind;
	double dist, min;

	for (i = 0; i < mw->rep.end; i++) {
		if (dominates(p->f, mw->rep.repo[i].f, mw->nobj)) {
			particle_copy(mw, &(mw->rep.repo[i]), p);
			return;
		}
		if (dominates(mw->rep.repo[i].f, p->f, mw->nobj)) {
			return;
		}
	}

	if (mw->rep.end < mw->nrepo) {
		particle_copy(mw, &(mw->rep.repo[mw->rep.end]), p);
		mw->rep.end++;
		return;
	}

	dist = repository_dist(mw, p);
	ind = 0;
	min = mw->rep.dist[0];
	for (i = 1; i < mw->rep.end; i++) {
		if (min > mw->rep.dist[i]) {
			min = mw->rep.dist[i];
			ind = i;
		}
	}
	if (dist > min) {
		particle_copy(mw, &(mw->rep.repo[ind]), p);
	}
}

void
repository_update(struct MOWAO *mw, struct Particle *p)
{
	int i, j;

	double obj[mw->nobj];
	if (mw->rep.end > 0) {
		for (i = 0; i < mw->nobj; i++)
			obj[i] = mw->rep.repo[0].f[i];
		mw->f(mw->rep.repo[0].dec, mw->ndec, mw->rep.repo[i].f,
		      mw->cb_arg);
		for (i = 0; i < mw->nobj; i++)
			if (obj[i] != mw->rep.repo[0].f[i])
				goto loop;
		goto next;
loop:
		for (i = 0; i < mw->rep.end; i++)
			mw->f(mw->rep.repo[i].dec, mw->ndec, mw->rep.repo[i].f,
			      mw->cb_arg);
	}
next:

	for (i = 0; i < mw->pop.end; i++) {
		for (j = i + 1; j < mw->pop.end; j++) {
			if (dominates(p[j].f, p[i].f, mw->nobj))
				goto nextp;
		}
		repository_add(mw, &p[i]);
nextp:;
	}

	for (i = 0; i < mw->pop.end; i++) {
		for (j = 0; j < mw->nobj; j++) {
			if (p[i].f[j] > mw->rep.max[j])
				mw->rep.max[j] = p[i].f[j];
		}
	}

	repository_check(mw);
}

double
repository_dist(struct MOWAO *mw, struct Particle *p)
{
	int i, j;
	double dist = INFINITY, d;

	for (i = 0; i < mw->rep.end; i++)
		mw->rep.dist[i] = INFINITY;

	for (i = 0; i < mw->rep.end; i++) {
		d = distance(p->f, mw->rep.repo[i].f, mw->nobj, mw->rep.max);
		if (d < dist)
			dist = d;
		for (j = i + 1; j < mw->rep.end; j++) {
			d = distance(mw->rep.repo[j].f, mw->rep.repo[i].f,
				     mw->nobj, mw->rep.max);
			if (d < mw->rep.dist[i]) {
				mw->rep.dist[i] = d;
			}
			if (d < mw->rep.dist[j])
				mw->rep.dist[j] = d;
		}
	}
	return dist;
}
