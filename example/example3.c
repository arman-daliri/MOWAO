#include <stdlib.h>
#include <math.h>

#include "../mowao/mowao.h"

/* ZDT2 */
void
f(double *x, int x_len, double *f, void *cb_arg)
{
	int i;
	double g, sum = 0;

	for (i = 1; i < x_len; i++)
		sum += x[i];
	f[0] = x[0];
	g = 1.0 + 9.0 * sum / (x_len - 1);
	f[1] = g * (1 - pow((f[0] * 1.0 / g), 2));
}

int
main(void)
{
	struct MOWAO *mw = mowao_new();
	int i;

	mw->maxiter = 100;
	mw->npop = 100;
	mw->nrepo = 80;
	mw->bond_radius = 0.001;
	mw->push = 0.5;
	mw->evaporate = 0.8;
	mw->coef = 1;

	mw->nobj = 2;
	mw->ndec = 10;

	mw->f = f;
	mowao_alloc(mw);

	for (i = 0; i < mw->ndec; i++) {
		mw->lb[i] = 0.0;
		mw->ub[i] = 1.0;
		mw->vlb[i] = -0.15;
		mw->vub[i] = 0.15;
	}

	mowao_init(mw);
	mowao_run(mw, 1);
	mowao_clean(mw);

	free(mw);
	mw = NULL;
	return 0;
}
