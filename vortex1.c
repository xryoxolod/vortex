#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <pthread.h>


// todo:
/* использование
 * vortex -h [H] -d ['[+|-](x1;y1) [+|-](x2;y2) ... координаты']
 * но вообще, формат
 * vortex -h [%f] -d ['[+|-] %f %f [+|-] %f %f ... координаты']
 */

// todo:
/* коды возврата
 * номер	-- описание
 */


/* настройки */
#define EF	2.07e-7		// Ф0	: Гс/см^2	// квант изменения магнитного поля

// размеры
#define X	5e-4		// ширина	: см
#define Y	5e-4		// высота	: см

// параметры образца
#define TC	84.0		// Tc	: K		// критическая температура
double tc = TC / 11606;		// Tc	: эВ		// критическая температура в Эв

#define L0	1.8e-5		// λ0	: см		// глубина проникновения при температуре T = 0 K
#define L	L0 * pow(1 - pow(t / tc, 3.3), -0.5)	// λ(T) : см	// глубина проникновения

#define K0	2e-7		// ξ0	: см		// длина когерентности при температуре T = 0 K
#define K	K0 * pow(1 - pow(t / tc, 3.3), -0.5)	// ξ(T) : см	// длинна когерентности

//#define DD	2.7e-8		// δ	: см		// высота слоя (условно по координите Z)
#define DD	2.7e-7		// δ	: см		// высота слоя (условно по координите Z)

#define U01	0.8		// U0	: эВ		// верхниий предел в распределении потенциалов дефектов
#define U02	0.1		// U0	: эВ		// нижний предел в распределении потенциалов дефектов


// todo: I
// распределение тока (если ток доступен)
//#define I	400		// I	: 


// распределение поля
//#define DISPX	2 * X / 9
//#define DISPY	2 * Y / 9
//#define MATHX	X / 2
//#define MATHY	Y / 2
//#define NORM	400.0
//#define PLUS	300.0
//#define H	exp(-((x - MATHX) * (x - MATHX) / (DISPX * DISPX) + (y - MATHY) * (y - MATHY) / (DISPY * DISPY)) / 2) / (2 * M_PI * DISPX * DISPY) / 12891550.390444 * NORM + PLUS

// постоянное поле
#define H	600		// H	: эрстед


// todo: I
// поле, связанное с током
double hi = 0;		// Hi	: эрстед


// распределение температуры
// Гаусс
//#define DISPX	2 * X / 9
//#define DISPY	2 * Y / 9
//#define MATHX	X / 2
//#define MATHY	Y / 2
//#define NORM	80.0 / 11606
//#define T	exp(-((x - MATHX) * (x - MATHX) / (DISPX * DISPX) + (y - MATHY) * (y - MATHY) / (DISPY * DISPY)) / 2) / (2 * M_PI * DISPX * DISPY) / 12891550.390444 * NORM

// Гаусс X
//#define DISPX	2 * X / 9
//#define MATHX	X / 2
//#define NORM	80.0 / 11606
//#define T	exp(-((x - MATHX) * (x - MATHX) / (DISPX * DISPX)) / 2) / (pow(2 * M_PI,0.5) * DISPX) / 3590.4805236128945580 * NORM

// Гаусс Y
//#define DISPY	2 * Y / 9
//#define MATHY	Y / 2
//#define NORM	80.0 / 11606
//#define T	exp(-((y - MATHY) * (y - MATHY) / (DISPY * DISPY)) / 2) / (pow(2 * M_PI,0.5) * DISPY) / 3590.4805236128945580 * NORM

// линейный потенциал
//#define T	(1 - x / X) * 40.0 / 11606	// T	: K

// постоянная температура
#define T	10.0 / 11606	// T	: K


// формулы расчёта энергий
double m1 = -EF / (4 * M_PI) * DD / (1.6 * 1e-12);
//double m11 = - 2 * M_PI * I / (1.6 * 1e-12 * 3 * 1e10); // todo: I
#define E1	m1 * h		// E1(H,λ,x,y)	: эВ	энергия от магнитного поля если все границы переодичны
#define E1X	m1 * (h * (1 - cosh((x - (X / 2)) / l) / cosh(X / (2 * l))) + hi * (sinh((x - (X / 2)) / l) / sinh(X / (2 * l)) + ((x > X / 2) ? -1 : 1) ))		// E1(H,λ,x,y)	: эВ	энергия от магнитного поля если границы не периодичны по X
#define E1Y	m1 * (h * (1 - cosh((y - (Y / 2)) / l) / cosh(Y / (2 * l))) + hi * (sinh((y - (Y / 2)) / l) / sinh(Y / (2 * l)) + ((y > Y / 2) ? -1 : 1) ))		// E1(H,λ,x,y)	: эВ	энергия от магнитного поля если границы не периодичны по Y
// todo: передалать энергию от магнитного поля
#define E1XY	m1 * h		// E1(H,λ,x)	: эВ	энергия от магнитного поля, если все границы не переодичны

double m2 = (EF * EF / (16 * M_PI * M_PI)) * DD / (1.6 * 1e-12);
#define E2	m2 * (log(l / xi) + 0.52) / (l * l)		// E2(λ,ξ)	: эВ	// собственная энергия вихря

double m3 = (EF * EF / (16 * M_PI * M_PI)) * DD / (1.6 * 1e-12);
#define E3	m3 * bessk(0, r / l) / (l * l)		// E3(r,λ) 	: эВ	// энергия взаимодействия типа вихрь-вихрь	: (194)

#define E4	u * exp(-r / (2 * xi)) / (r / xi + 1)		// E4(r,ξ)	: эВ	// энергия взаимодействия типа примесь-вихрь


// примеси
#define NP	0		// количество примесей	: штук
#define RP	0.04		// отвечает за хаотичность примесей	: безразм. (чем коэфициент меньше, тем упорядоченней будут вести себя вихри)
#define CP	1		// 0 - считать взаимодействие всех примесей с вихрем; 1 - считать взаимодействие ТОЛЬКО ближайшей примеси с вихрем


// расстояния
#define DP	0.1		// максимальный радиус перемещения вихря (от λ)	: см
#define DS	10 * K0		// максимальное расстояние на котором допустима аннигиляция	: см
#define DR	K0		// минимальное расстояние на которое подходят вихри	: см
#define BRN	L0		// максимальное расстояние от стенки, на котором рождаются вихри	: см
#define DRN	0.4 * L0	// минимальное расстояние от стенки, на котором рождаются вихри	: см

#define TA	1		// τ (коэффициент для подкручивания вероятности)


// шаги
// todo: если D != -1, то криво работает
#define D	-1		// длина времени измерений (сколько програмных шагов в одном физическом); при -1 будет пытаться ориентироваться на количество вихрей (ЭТО НЕ БУДЕТ ФИКСИРОВАННАЯ ВЕЛИЧИНА)
//#define DM	4		// минимальное колличество програмных шагов в одном физическом (используется если D == -1)
#define DM	4000		// минимальное колличество програмных шагов в одном физическом (используется если D == -1)
#define ST	8		// максимальное колличество программных шагов, которое выполяет ядро (при -1 неограничено) // fail: не работает
#define RND	4		// коэффициент, определяющий вероятность выбора действия в програмном шаге (добавить, убрать, передвинуть)


// fail: не работает
/*
#define EE	1		// определяет то, как считать стационарные энергии (E1 + E2 + E4), если 0, то будет БЫСТРО брать ПРИБЛИЖЁННЫЕ значения из массива, если 1, то будет МЕДЛЕННО высчитывать ТОЧНЫЕ значения в течение всего расчсёт
unsigned long long ex = (int) (X / (K0 / 4)),	// количество разбиений системы по направлению X при приближённом расчёте энергиий
	      ey = (int)(Y / (K0 / 4));		// количество разбиений системы по направлению Y при приближённом расчёте энергиий
*/


// границы
#define XG	1		// определяет характер границ
// 0	-- все границы периодеческие
// 1	-- периодичны границы по Y
// 2	-- периодичны границы по X
// 3	-- нет периодичных границ // fail: не работает так как надо

// выход из программы, погрешность
#define G	0.5		// часть массива с конца, берётся для разбиения значений в погрешности
#define P	0.05		// допустимая погрешность (0.05 ~ 0.1) (лучше ставить 0.1)
/* --- */


// todo:
/* переменные
 * h	-- H	// внешнее магнитное поле	: эрстеды
 * t	-- T	// температура			: эВ
 * r	-- r	// расстояние			: см
 * x	-- x	// координата по х		: см
 * y	-- y	// координата по у		: см
 * l	-- λ	// глубина проникновения	: см
 * xi	-- ξ	// длинна когерентности		: см
 */


// функции бесселя из bessel.c
extern double bessj(int,double);
extern double bessy(int,double);
extern double bessi(int,double);
extern double bessk(int,double);


struct vortex
{
	double x;	// координата по x
	double y;	// координата по y
	char zn;	// направление вихря (o -- на нас ('+' в энергиях) / x -- в нас ('-' в энергиях))
	double t;	// T, температура в точке
	double h;	// H, поле в точке
	double l;	// λ, глубина проникновения в точке
	double xi;	// ξ, длина когерентности в точке
	double e;	// E1 + E2 + E4, все стационарные энергии
} **v;			// массив, содержащий вихри всего расчёта
struct prim
{
	double x;	// координата по x
	double y;	// координата по y
	double u;	// потенциал дефекта
} *pr;			// массив, содержащий примеси
int *nn;		// массив, содержащий колличество вихрей от физического шага (nn[номер шага] = колличество вихрей)
double *e,		// массив, содержащий полную энергию системы от физического шага (e[номер шага] = энергия)      : эВ
       e0,		// изменение энергии происходящее на каждом програмном шаге
       *me,		// массив, со средними значениями энергий с последних G * nm шагов
       **ee,		// двумерный массив, хранит примерные значения стационарных энергий (e[x][y] = примерное значение энергии)
       pog;		// текущая погрешность
unsigned long long k,	// програмные шаги
	      dk = 0,	// количество програмных шагов до записи следующего физического шага (счётчик)
	      nm,	// физические шаги
	      d;	// номер следующего физического шага в програмном понимании



// todo:
/* обработчик ошибок */
void erro(char *op, ...) {
	va_list f;
	va_start(f, op);

	int i = va_arg(f, int);
	switch (i)
	{
		case 1:
		case 2:
		case 3:
		case 4:
			printf("ОШИБКА ВВОДА: ");
			if (1 < i && i < 4) printf("КООРДИНАТЫ>> ");
			printf("%s\n\tдля просмотра справки вызовите программу с ключом -help.\n", op);
			break;
		default:
			printf("%s\n", op);
	}

	va_end(f);
	exit(i);
}

/* модуль */
double m(double i) {
        return i > 0 ? i : -i;
}

/* меняет местами два double */
void swapch(double *a, double *b) {
	double c = *a; *a = *b; *b = c;
}

// todo: передалать поле
int qq = 0, tm = 0;
/* магнитное поле H(x,y) */
double pole(double x, double y) {
	return (!tm) ? H : H;
}

// todo: передалать температуру
/* температура T(x,y) */
double temp(double x, double y) {
	return (!tm) ? T : 5.0 / 11606;
}

/* глубина проникновения λ(T) */
double lamb(double t) {
	return L;
}

/* длина когерентности ξ(T) */
double ksii(double t) {
	return K;
}

/* энергия взаимодействия с магнитным полем E1(H,λ,x,y) */
double e1(double h, double l, double x, double y) {
#if XG == 0
	return E1;
#endif
#if XG == 1
	return E1Y;
#endif
#if XG == 2
	return E1X;
#endif
#if XG == 3
	return E1XY;
#endif
}

/* собственая энергия вихря E2(λ,ξ) */
double e2(double l, double xi) {
	return E2;
}

/* энергия взаимодействия вихря с другим вихрем E3(r,λ) */
double e3(double r, double l) {
	return E3;
}

/* энергия взаимодействия вихря с примесью E4(r,ξ,U) */
double e4(double r, double xi, double u) {
	return E4;
}

/* определяет минимальное расстояние между двумя точками с учётом стенок и границ системы */
double rrr(double x1, double y1, double x2, double y2) {
	double xi, yi, x, y;

	xi = m(x1 - x2); yi = m(y1 - y2);
#if XG & 2
	x = xi;
#else
	x = (xi < X - xi) ? xi : X - xi;
#endif
#if XG & 1
	y = yi;
#else
	y = (yi < Y - yi) ? yi : Y - yi;
#endif

	return pow(x * x + y * y, 0.5);
}

/* определяет минимальное расстояние между двумя точками без учёта стенок и границ системы */
double rrb(double x1, double y1, double x2, double y2) {
	double x = m(x1 - x2), y = m(y1 - y2);

	return pow(x * x + y * y, 0.5);
}

/* проверяет, есть ли в vortex *a вихрь рядом с координатами (x,y) (vortex *a имеет размер n) */
char check(struct vortex *a, int n, double x, double y) {
	for (int i = 0; i < n; i++)
		if (rrr(a[i].x,a[i].y,x,y) <= DR) return 1;

	return 0;
}

/* меняет местами два вихря */
void swapvi(struct vortex *a, int i, int j) {
	double s;

	s = a[i].x; a[i].x = a[j].x; a[j].x = s;
	s = a[i].y; a[i].y = a[j].y; a[j].y = s;
	s = a[i].zn; a[i].zn = a[j].zn; a[j].zn = s;
	s = a[i].t; a[i].t = a[j].t; a[j].t = s;
	s = a[i].h; a[i].h = a[j].h; a[j].h = s;
	s = a[i].l; a[i].l = a[j].l; a[j].l = s;
	s = a[i].xi; a[i].xi = a[j].xi; a[j].xi = s;
	s = a[i].e; a[i].e = a[j].e; a[j].e = s;
}

/* копирует struct vortex *b в память и возвращает ссылку на копию (struct vortex *b имеет размер n) */
struct vortex *cpo(struct vortex *b, int n) {
	struct vortex *a;

	a = malloc(n * sizeof(struct vortex));
	for (int i = 0; i < i; i++)
	{
		a[i].x = b[i].x;
		a[i].y = b[i].y;
		a[i].zn = b[i].zn;
		a[i].t = b[i].t;
		a[i].h = b[i].h;
		a[i].l = b[i].l;
		a[i].xi = a[i].xi;
		a[i].e = b[i].e;
	}

	return a;
}


/* --- */


/* считает энергию от магнитного поля всей vortex *a (vortex *a имеет размер n) */
double tote1(struct vortex *a, int n) {
	double v = 0, f = 0, s, t;

	for (int i = 0; i < n; i++)
	{
		s = ((a[i].zn == 'o') ? 1 : -1) * e1(a[i].h,a[i].l,a[i].x,a[i].y) - f;
		t = v + s;
		f = (t - v) - s;
		v = t;
	}

	return v;
}

/* считает внутреннию энергию всех вихрей всей vortex *a (vortex *a имеет размер n) */
double tote2(struct vortex *a, int n) {
	double v = 0, f = 0, s, t;

	for (int i = 0; i < n; i++)
	{
		s = e2(a[i].l,a[i].xi) - f;
		t = v + s;
		f = (t - v) - s;
		v = t;
	}

	return v;
}

/* считает энергию взаимодействий типа вихрь-вихрь последней в vortex *a частицы (vortex *a имеет размер n) */
double eone3(struct vortex *a, int n) {
	double v = 0, f = 0, s, t, r;

	for (int i = 0; i < n - 1; i++)
	{
		r = rrr(a[n - 1].x,a[n - 1].y,a[i].x,a[i].y);
		if (a[i].zn == a[n - 1].zn)
			s = e3(r,a[n - 1].l) + e3(r,a[i].l) - f;
		else
			s = -e3(r,a[n - 1].l) - e3(r,a[i].l) - f;
		t = v + s;
		f = (t - v) - s;
		v = t;

#if XG & 1
		// todo:
		// if (a[i].x < X / 2)

		r = rrb(a[n - 1].x,a[n - 1].y,-a[i].x,a[i].y);
		if (a[i].zn == a[n - 1].zn)
			s = e3(r,a[n - 1].l) + e3(r,a[i].l) - f;
		else
			s = -e3(r,a[n - 1].l) - e3(r,a[i].l) - f;
		t = v + s;
		f = (t - v) - s;
		v = t;

		r = rrb(a[n - 1].x,a[n - 1].y,2 * X - a[i].x,a[i].y);
		if (a[i].zn == a[n - 1].zn)
			s = e3(r,a[n - 1].l) + e3(r,a[i].l) - f;
		else
			s = -e3(r,a[n - 1].l) - e3(r,a[i].l) - f;
		t = v + s;
		f = (t - v) - s;
		v = t;
#endif
#if XG & 2
		// todo:
		//if (a[i].y < Y / 2)

		r = rrb(a[n - 1].x,a[n - 1].y,a[i].x,-a[i].y);
		if (a[i].zn == a[n - 1].zn)
			s = e3(r,a[n - 1].l) + e3(r,a[i].l) - f;
		else
			s = -e3(r,a[n - 1].l) - e3(r,a[i].l) - f;
		t = v + s;
		f = (t - v) - s;
		v = t;

		r = rrb(a[n - 1].x,a[n - 1].y,a[i].x,2 * Y - a[i].y);
		if (a[i].zn == a[n - 1].zn)
			s = e3(r,a[n - 1].l) + e3(r,a[i].l) - f;
		else
			s = -e3(r,a[n - 1].l) - e3(r,a[i].l) - f;
		t = v + s;
		f = (t - v) - s;
		v = t;
#endif
	}
#if XG & 1
	// todo:
	// if (a[i].x < X / 2)

	r = rrb(a[n - 1].x,a[n - 1].y,-a[n - 1].x,a[n - 1].y);
	s = -e3(r,a[n - 1].l) - f;
	t = v + s;
	f = (t - v) - s;
	v = t;

	r = rrb(a[n - 1].x,a[n - 1].y,2 * X - a[n - 1].x,a[n - 1].y);
	s = -e3(r,a[n - 1].l) - f;
	t = v + s;
	f = (t - v) - s;
	v = t;
#endif
#if XG & 2
	// todo:
	//if (a[i].y < Y / 2)

	r = rrb(a[n - 1].x,a[n - 1].y,a[n - 1].x,-a[n - 1].y);
	s = -e3(r,a[n - 1].l) - f;
	t = v + s;
	f = (t - v) - s;
	v = t;

	r = rrb(a[n - 1].x,a[n - 1].y,a[n - 1].x,2 * Y - a[n - 1].y);
	s = -e3(r,a[n - 1].l) - f;
	t = v + s;
	f = (t - v) - s;
	v = t;
#endif

	return v;
}

/* считает всю энергию взаимодействий типа вихрь-вихрь в vortex *a (vortex *a имеет размер n) */
double tote3(struct vortex *a, int n) {
	double r, v = 0, f = 0, s, t;
	int i, j;

	for (j = 0; j < n; j++)
	{
		for (i = 0; i < j; i++)
		{
			r = rrr(a[j].x,a[j].y,a[i].x,a[i].y);
			if (a[i].zn == a[j].zn)
				s = e3(r,a[j].l) + e3(r,a[i].l) - f;
			else
				s = -e3(r,a[j].l) - e3(r,a[i].l) - f;
			t = v + s;
			f = (t - v) - s;
			v = t;

#if XG & 1
			// todo:
			// if (a[i].x < X / 2)

			r = rrb(a[j].x,a[j].y,-a[i].x,a[i].y);
			if (a[j].zn == a[i].zn)
				s = e3(r,a[j].l) + e3(r,a[i].l) - f;
			else
				s = -e3(r,a[j].l) - e3(r,a[i].l) - f;
			t = v + s;
			f = (t - v) - s;
			v = t;

			r = rrb(a[j].x,a[j].y,2 * X - a[i].x,a[i].y);
			if (a[i].zn == a[j].zn)
				s = e3(r,a[j].l) + e3(r,a[i].l) - f;
			else
				s = -e3(r,a[j].l) - e3(r,a[i].l) - f;
			t = v + s;
			f = (t - v) - s;
			v = t;
#endif
#if XG & 2
			// todo:
			//if (a[i].y < Y / 2)

			r = rrb(a[j].x,a[j].y,a[i].x,-a[i].y);
			if (a[j].zn == a[i].zn)
				s = e3(r,a[j].l) + e3(r,a[i].l) - f;
			else
				s = -e3(r,a[j].l) - e3(r,a[i].l) - f;
			t = v + s;
			f = (t - v) - s;
			v = t;

			r = rrb(a[j].x,a[j].y,a[i].x,2 * Y - a[i].y);
			if (a[i].zn == a[j].zn)
				s = e3(r,a[j].l) + e3(r,a[i].l) - f;
			else
				s = -e3(r,a[j].l) - e3(r,a[i].l) - f;
			t = v + s;
			f = (t - v) - s;
			v = t;
#endif
		}
#if XG & 1
		// todo:
		// if (a[i].x < X / 2)

		r = rrb(a[j].x,a[j].y,-a[j].x,a[j].y);
		s = -e3(r,a[j].l) - f;
		t = v + s;
		f = (t - v) - s;
		v = t;

		r = rrb(a[j].x,a[j].y,2 * X - a[j].x,a[j].y);
		s = -e3(r,a[j].l) - f;
		t = v + s;
		f = (t - v) - s;
		v = t;
#endif
#if XG & 2
		// todo:
		//if (a[i].y < Y / 2)

		r = rrb(a[j].x,a[j].y,a[j].x,-a[j].y);
		s = -e3(r,a[j].l) - f;
		t = v + s;
		f = (t - v) - s;
		v = t;

		r = rrb(a[j].x,a[j].y,a[j].x,2 * Y - a[j].y);
		s = -e3(r,a[j].l) - f;
		t = v + s;
		f = (t - v) - s;
		v = t;
#endif
	}

	return v;
}

/* считает энергию взаимодействий типа вихрь-примесь n-го вихря в vortex *a */
double eone4(struct vortex *a, int n) {
	double r, v = 0, f = 0, s, t;

	for (int i = 0; i < NP; i++)
	{
		s = e4(rrr(a[n].x,a[n].y,pr[i].x,pr[i].y),a[n].xi,pr[i].u) - f;
		t = v + s;
		f = (t - v) - s;
		v = t;
	}

	return v;
}

/* считает всю энергию взаимодействий типа вихрь-примесь в vortex *a (vortex *a имеет размер n) */
double tote4(struct vortex *a, int n) {
	double r, v = 0, f = 0, s, t;
	int i, j;

	for (i = 0; i < NP; i++)
		for (j = 0; j < n; j++)
		{
			s = e4(rrr(a[j].x,a[j].y,pr[i].x,pr[i].y),a[j].xi,pr[i].u) - f;
			t = v + s;
			f = (t - v) - s;
			v = t;
		}

	return v;
}

/* считает все стационарные энергии, записаные заранее в vortex *a (vortex *a имеет размер n) */
double tote124(struct vortex *a, int n) {
	double r, v = 0, f = 0, s, t;
	int i, j;

	for (i = 0; i < n; i++)
	{
		s = a[i].e - f;
		t = v + s;
		f = (t - v) - s;
		v = t;
	}

	return v;
}


/* --- */


/* программа показывает статус при нажатии 's' */
void show() {
	system("clear");
	printf("шагов: %d / %llu\nвихрей: %d\ntotal: %.16f\n+++++: %.16f\nпогрешность: %5.2f / %5.2f %\n", nm, k, nn[nm], tote1(v[nm],nn[nm]) + tote2(v[nm],nn[nm]) + tote3(v[nm],nn[nm]) + tote4(v[nm],nn[nm]), e[nm], 100 * pog, 100 * (double)P);
	fflush(stdout);
}

/* создаёт случайный новый вихрь в vortex *a (vortex *a имеет размеры n) */
struct vortex *new(struct vortex *a, int *n) {
	char zn = rand() % 2 ? 'o' : 'x';
	double x, y;

#if XG & 2
	x = (double)rand() / RAND_MAX * (BRN - DRN) + (rand() % 2 ? DRN : X - L0);
#else
	x = (double)rand() / RAND_MAX * X;
#endif
#if XG & 1
	y = (double)rand() / RAND_MAX * (BRN - DRN) + (rand() % 2 ? DRN : Y - L0);
#else
	y = (double)rand() / RAND_MAX * Y;
#endif

	if (check(a,*n,x,y)) return 0;

	a = realloc(a, (*n + 1) * sizeof(struct vortex));

	a[*n].x = x;
	a[*n].y = y;
	a[*n].zn = zn;
	a[*n].t = temp(x,y);
	a[*n].h = pole(x,y);
	a[*n].l = lamb(a[*n].t);
	a[*n].xi = ksii(a[*n].t);
	a[*n].e = ((a[*n].zn == 'o') ? 1 : -1) * e1(a[*n].h,a[*n].l,a[*n].x,a[*n].y) + e2(a[*n].l,a[*n].xi) + eone4(a,*n);

	*n = *n + 1;

	return a;
}


char qqq = 0;
double ff = 0, fs = 0;
/* условие остановки */
char stopr(double ee) {
	int i;
	double f, s, g;


	s = ee - ff;
	g = e[nm] + s;
	f = (g - e[nm]) - s;
	e[nm] = g;

#if D < 0
	if (k > d + dk)
	{
		dk = k;
		d = (nn[nm] < DM) ? DM : nn[nm];
#else
	if (k > D + dk)
	{
#endif
		nm++;
		v = realloc(v, (nm + 1) * sizeof(struct vortex*));
		v[nm] = malloc(nn[nm - 1] * sizeof(struct vortex));
		nn = realloc(nn, (nm + 1) * sizeof(int));
		nn[nm] = nn[nm - 1];
		for (i = 0; i < nn[nm]; i++)
		{
			v[nm][i].x = v[nm - 1][i].x;
			v[nm][i].y = v[nm - 1][i].y;
			v[nm][i].zn = v[nm - 1][i].zn;
			v[nm][i].t = v[nm - 1][i].t;
			v[nm][i].h = v[nm - 1][i].h;
			v[nm][i].l = v[nm - 1][i].l;
			v[nm][i].xi = v[nm - 1][i].xi;
			v[nm][i].e = v[nm - 1][i].e;
		}
		e = realloc(e, (nm + 1) * sizeof(double));
		e[nm] = e[nm - 1];

		me = realloc(me, (nm + 1) * sizeof(double));
		for (me[nm - 1] = i = 0; i < nm; i++)
		{
			s = e[i] - fs;
			g = me[nm - 1] + s;
			fs = (g - me[nm - 1]) - s;
			me[nm - 1] = g;
		}
		me[nm - 1] = me[nm - 1] / nm;

		// s -- минимальное значение
		// g -- максимальное значение
		for (s = g = me[nm - 1], i = (int)((1 - G) * nm); i < nm; i++)
		{
			if (s > me[i]) s = me[i];
			if (g < me[i]) g = me[i];
		}
		pog = m(((g - s) / 2) / e[nm]);
		if (!pog) pog = 1;

		// условия выхода
//		if (pog <= P && !q)
//		q = 2 * nm;
//		if (q && q < nm)
//			if (pog > P)
//				q = 0;
//			else
//				return 1;

//		if (nm > 60)
//			qqq = 1;

		if (pog <= P && !qq && !tm)
			qq = 2 * nm;
		if (qq && qq < nm && !tm)
			if (pog > P)
				qq = 0;
			else {tm = 1; qq = 0; show();}
		if (!qq && tm)
		{
			qq = 2 * nm;
			for (i = 0; i < nn[nm]; i++)
			{
				v[nm][i].t = temp(v[nm][i].x,v[nm][i].y);
				v[nm][i].h = pole(v[nm][i].x,v[nm][i].y);
				v[nm][i].l = lamb(v[nm][i].t);
				v[nm][i].xi = ksii(v[nm][i].t);
				v[nm][i].e = ((v[nm][i].zn == 'o') ? 1 : -1) * e1(v[nm][i].h,v[nm][i].l,v[nm][i].x,v[nm][i].y) + e2(v[nm][i].l,v[nm][i].xi) + eone4(v[nm],i);
			}
			e[nm] = tote1(v[nm],nn[nm]) + tote2(v[nm],nn[nm]) + tote3(v[nm],nn[nm]) + tote4(v[nm],nn[nm]);
		}
		if (qq && qq < nm && tm)
			qqq = 1; // заканчиваем работу

		if (nm > 10000)
			qqq = 1; // слишком много
		// todo: иногда погрешность может быть больше 100
//		if (isnan(pog) || pog > 100)
//			qqq = 1; // заканчиваем работу
	}

	return 0;
}

/* ветка основного алгоритма */
void *three() {
	int i;
	double g, x, y, t, h, l, xi;


	e0 = 0;
	k++;

	i = rand() % RND;
	if (i == 1 || nn[nm] != 0)
	{
		switch (i)
		{
			case 0: // удаление
				i = rand() % nn[nm];

#if XG & 2
	if ((v[nm][i].y <= Y - L0) && (v[nm][i].y >= L0)) break;
#endif
#if XG & 1
	if ((v[nm][i].x <= X - L0) && (v[nm][i].x >= L0)) break;
#endif

				swapvi(v[nm],i,nn[nm] - 1);

				e0 = -(eone3(v[nm],nn[nm]) + v[nm][nn[nm] - 1].e);

				g = nn[nm] / TA * (double)nn[nm] * exp(-e0 / v[nm][nn[nm] - 1].t);
				if (g < 1)
				{
					if ((double)rand() / RAND_MAX <= g)
					{
						v[nm] = realloc(v[nm], nn[nm] * sizeof(struct vortex));
						nn[nm]--;
					}
					else e0 = 0;
				}
				else
				{
					v[nm] = realloc(v[nm], nn[nm] * sizeof(struct vortex));
					nn[nm]--;
				}
				break;
			case 1: // создание
				struct vortex *a;

//				v[nm] = new(v[nm],&nn[nm]);
				a = new(v[nm],&nn[nm]);
				if (!a) break;
				v[nm] = a;

				e0 = eone3(v[nm],nn[nm]) + v[nm][nn[nm] - 1].e;

				g = (double)TA / (nn[nm] + 1) * exp(-e0 / v[nm][nn[nm] - 1].t);
				if (g < 1)
				{
					if ((double)rand() / RAND_MAX > g)
					{
						v[nm] = realloc(v[nm], nn[nm] * sizeof(struct vortex));
						nn[nm]--;
						e0 = 0;
					}
				}
				break;
			case 2: // анигиляция
				// fail: вихри не анигилируют
				int j, xx, yy;

				for (i = xx = yy = 0; i < nn[nm]; i++)
					if (v[nm][i].zn == 'o') yy++;
					else xx++;

				if (xx && yy)
				{
					i = rand() % xx;
					j = rand() % yy;
					for (int tt = yy = xx = 0; tt < nn[nm]; tt++)
						if (v[nm][tt].zn == 'o')
						{
							if (yy == j) j = tt;
							yy++;
						}
						else
						{
							if (xx == i) i = tt;
							xx++;
						}

					if ("возможна аннигиляция одного вихря", 421) // erro 421

					if (rrr(v[nm][i].x,v[nm][i].y,v[nm][j].x,v[nm][j].y) <= DS)
					{
						e0 = -(v[nm][i].e + v[nm][j].e);

						swapvi(v[nm],i,nn[nm] - 1);
						e0 = e0 - eone3(v[nm],nn[nm]);

						v[nm] = realloc(v[nm], nn[nm] * sizeof(struct vortex));
						nn[nm]--;

						swapvi(v[nm],(j == nn[nm]) ? i : j,nn[nm] - 1);
						e0 = e0 - eone3(v[nm],nn[nm]);

						v[nm] = realloc(v[nm], nn[nm] * sizeof(struct vortex));
						nn[nm]--;
					}
				}
				break;
			default: // перемещение
				i = rand() % nn[nm];

				do {
					x = (double)rand() / RAND_MAX * 2 * DP * v[nm][i].l - DP * v[nm][i].l;
					y = (double)rand() / RAND_MAX * 2 * DP * v[nm][i].l - DP * v[nm][i].l;
				} while (x * x + y * y > DP * v[nm][i].l);
#if XG & 2
				if (v[nm][i].x + x < 0)
					break;
				else
					if (v[nm][i].x + x >= X)
						break;
					else
						x = v[nm][i].x + x;
#else
				x = (v[nm][i].x + x < 0) ? (X + v[nm][i].x + x) : ((v[nm][i].x + x >= X) ? (v[nm][i].x + x - X) : (v[nm][i].x + x));
#endif
#if XG & 1
				if (v[nm][i].y + y < 0)
					break;
				else
					if (v[nm][i].y + y >= Y)
						break;
					else
						y = v[nm][i].y + y;
#else
				y = (v[nm][i].y + y < 0) ? (Y + v[nm][i].y + y) : ((v[nm][i].y + y >= Y) ? (v[nm][i].y + y - Y) : (v[nm][i].y + y));
#endif
				if (check(v[nm],nn[nm],x,y)) break;

				swapvi(v[nm],i,nn[nm] - 1);

				e0 = eone3(v[nm],nn[nm]) + v[nm][nn[nm] - 1].e;
				t = v[nm][nn[nm] - 1].t; v[nm][nn[nm] - 1].t = temp(x,y);
				h = v[nm][nn[nm] - 1].h; v[nm][nn[nm] - 1].h = pole(x,y);
				g = v[nm][nn[nm] - 1].x; v[nm][nn[nm] - 1].x = x; x = g;
				g = v[nm][nn[nm] - 1].y; v[nm][nn[nm] - 1].y = y; y = g;
				l = v[nm][nn[nm] - 1].l; v[nm][nn[nm] - 1].l = lamb(v[nm][nn[nm] - 1].t);
				xi = v[nm][nn[nm] - 1].xi; v[nm][nn[nm] - 1].xi = ksii(v[nm][nn[nm] - 1].t);
				g = v[nm][nn[nm] - 1].e; v[nm][nn[nm] - 1].e = ((v[nm][nn[nm] - 1].zn == 'o') ? 1 : -1) * e1(v[nm][nn[nm] - 1].h,v[nm][nn[nm] - 1].l,v[nm][nn[nm] - 1].x,v[nm][nn[nm] - 1].y) + e2(v[nm][nn[nm] - 1].l,v[nm][nn[nm] - 1].xi) + eone4(v[nm],nn[nm] - 1);
				e0 = eone3(v[nm],nn[nm]) + v[nm][nn[nm] - 1].e - e0;
				if (e0 > 0)
				{
					if ((double)rand() / RAND_MAX > (double)exp(-e0 / v[nm][nn[nm] - 1].t))
					{
						e0 = 0;
						v[nm][nn[nm] - 1].x = x;
						v[nm][nn[nm] - 1].y = y;
						v[nm][nn[nm] - 1].t = t;
						v[nm][nn[nm] - 1].h = h;
						v[nm][nn[nm] - 1].l = l;
						v[nm][nn[nm] - 1].xi = xi;
						v[nm][nn[nm] - 1].e = g;
					}
				}
		}
	}
	stopr(e0);
	pthread_exit(0);
}


/* --- */


/* распределяет в системе примеси */
void *defect() {
	int i;
	double w = 1.32471795724474602596, a1, a2, xi, yi;

	a1 = 1.0 / w;
	a2 = 1.0 / (w * w);

	pr = malloc(NP * sizeof(struct prim));
	for (i = 0; i < NP; i++)
	{
		pr[i].x = ((double)rand() - (double)rand()) / RAND_MAX * RP;
		pr[i].y = ((double)rand() - (double)rand()) / RAND_MAX * RP;
		w = ((0.5 + a1 * i) - (int)(0.5 + a1 * i));
		pr[i].x = ((pr[i].x + w > 1) ? (pr[i].x + w - 1) : ((pr[i].x + w < 0) ? (pr[i].x + w + 1) : (pr[i].x + w))) * X;
		w = ((0.5 + a2 * i) - (int)(0.5 + a2 * i));
		pr[i].y = ((pr[i].y + w > 1) ? (pr[i].y + w - 1) : ((pr[i].y + w < 0) ? (pr[i].y + w + 1) : (pr[i].y + w))) * X;
		pr[i].u = (double)rand() / RAND_MAX * m(U01 - U02) + U02;
	}
}


/* принимает значения по ключам при вызове программы */
void take(int argc, char *argv[]) {
	double f, s, t;
	int i, j, x, y, p, q;

	v = malloc(sizeof(struct vortex *));
	v[0] = malloc(sizeof(struct vortex));
	nn = malloc(sizeof(int));
	v[0][0].zn = 'o';

	for (i = 1; i < argc; i++)
	{
		// todo: help
		if ((argv[i][0] == '-' && argv[i][1] == 'h' && argv[i][2] == 'e' && argv[i][3] == 'l' && argv[i][4] == 'p' && argv[i][5] == '\0') || (argv[i][0] == '-' && argv[i][1] == '-' && argv[i][2] == 'h' && argv[i][3] == 'e' && argv[i][4] == 'l' && argv[i][5] == 'p' && argv[i][6] == '\0') || (argv[i][0] == 'h' && argv[i][1] == 'e' && argv[i][2] == 'l' && argv[i][3] == 'p' && argv[i][4] == '\0')) // help
		{
			printf("используются следующие ключи:\
				\n  -d    -- для задания начальной конфигурации расположения вихрей\
				\n  -p    -- для задания начальной конфигурации расположения примесей\
				\n  -help -- для выхова этой справки\
				\n");

			erro("запрошен help", 0); // erro 0
		}
		if (argv[i][0] == '-' && argv[i][1] == 'd' && argv[i][2] == '\0') // вихри
		{
			i++; y = x = p = j = 0;
			while (argv[i][j] != '\0' && argv[i][j] != '|')
			{
				if (argv[i][j] == '-')
					if (y == 0 || y == 2)
					{
						x = 1;
						j++;
						continue;
					}
					else erro("минус в неожиданом месте", 3); // erro 3
				if (argv[i][j] == '.')
					if (!p) {p = 1; j++; continue;}
					else erro("слишком много точек во вводе", 1); // erro 1
				if (argv[i][j] == 'x') {v[nm][nn[nm]].zn = 'x'; j++; continue;}
				if ('0' <= argv[i][j] && '9' >= argv[i][j])
				{
					if (y == 0) y = 1;
					if (y == 1)
						if (p)
						{
							v[nm][nn[nm]].x = v[nm][nn[nm]].x + (argv[i][j] - '0') * pow(10,-p);
							p++;
						}
						else
							v[nm][nn[nm]].x = 10 * v[nm][nn[nm]].x + argv[i][j] - '0';
					if (y == 2) y = 3;
					if (y == 3)
						if (p)
						{
							v[nm][nn[nm]].y = v[nm][nn[nm]].y + (argv[i][j] - '0') * pow(10,-p);
							p++;
						}
						else
							v[nm][nn[nm]].y = 10 * v[nm][nn[nm]].y + argv[i][j] - '0';
				}
				else
				{
					if (y == 3)
					{
						if (x) v[nm][nn[nm]].y = -v[nm][nn[nm]].y;
						p = x = y = 0;

						nn[nm]++;
						v[nm] = realloc(v[nm], (nn[nm] + 1) * sizeof(struct vortex));
						v[nm][nn[nm]].x = 0;
						v[nm][nn[nm]].y = 0;
						v[nm][nn[nm]].zn = 'o';
					}
					if (y == 1)
					{
						y = 2;
						if (x) v[nm][nn[nm]].x = -v[nm][nn[nm]].x;
						p = x = 0;
					}
				}
				j++;
			}
		}
		// todo:
		// if (argv[i][0] == '-' && argv[i][1] == 'p' && argv[i][2] == '\0') // примеси
		// {
		// }
	}
	nm++;
	v = realloc(v, (nm + 1) * sizeof(struct vortex*));
	v[nm] = malloc(nn[nm - 1] * sizeof(struct vortex));
	nn = realloc(nn, (nm + 1) * sizeof(int));
	nn[nm] = nn[nm - 1];
	for (i = 0; i < nn[nm]; i++)
	{
		v[nm][i].x = v[nm - 1][i].x;
		v[nm][i].y = v[nm - 1][i].y;
		v[nm][i].zn = v[nm - 1][i].zn;
		v[nm][i].t = v[nm - 1][i].t = temp(v[nm][i].x,v[nm][i].y);
		v[nm][i].h = v[nm - 1][i].h = pole(v[nm][i].x,v[nm][i].y);
		v[nm][i].l = v[nm - 1][i].l = lamb(v[nm][i].t);
		v[nm][i].xi = v[nm - 1][i].xi = ksii(v[nm][i].t);
		v[nm][i].e = v[nm - 1][i].e = ((v[nm][i].zn == 'o') ? 1 : -1) * e1(v[nm][i].h,v[nm][i].l,v[nm][i].x,v[nm][i].y) + e2(v[nm][i].l,v[nm][i].xi) + eone4(v[nm],i);
	}
}

/* сохраняет всю информацию на диск в файл с названием *ss */
void save(char *ss) {
	int i, j;
	FILE *out;

	out = fopen(ss, "w");
	fprintf(out, "δE = %.4f эВ\n* : %.16f\n", m(pog * e[nm]), e[0]);
	for (i = 0; i < nn[0]; i++)
		fprintf(out, "%c(%.8f;%.8f) ", v[0][i].zn, v[0][i].x, v[0][i].y);
	fprintf(out, "\n");
	for (i = 0; i < nm; i++)
	{
		fprintf(out, "%lld : %.16f\n", (unsigned long long)i, e[i]);
		for (j = 0; j < nn[i]; j++)
			fprintf(out, "%c(%.8f;%.8f) ", v[i][j].zn, v[i][j].x, v[i][j].y);
		fprintf(out, "\n");
	}
	fflush(out);
	fclose(out);
}


/* --- */


/* |vvv| */
/* печатает файл с местоположением примесей */
void *primes() {
	FILE *rpim;

	rpim = fopen("prim.txt", "w");
	for (int i = 0; i < NP; i++)
		fprintf(rpim, "%.16f;%.16f\n", pr[i].x, pr[i].y);

	fflush(rpim);
	fclose(rpim);
}

// todo: подумать над разбиением
/* печатает файл с распределением температуры (T) в системе */
void *pritemp() {
	unsigned long long ex = 2500,	// количество разбиений системы по направлению X
		      ey = 2500;	// количество разбиений системы по направлению Y
	int i, j;
	double x = X / (2 * ex), dx = X / ex, y = Y / (2 * ey), dy = Y / ey;
	FILE *etem;

	etem = fopen("T.txt", "w");
	for (i = 0; i <= ex; i++)
	{
		for (j = 0; j <= ey; j++)
			fprintf(etem, "%.12f ", temp(x + dx * i,y + dy * j));
		fprintf(etem, "\n");
	}
	fflush(etem);
	fclose(etem);
}

// todo: подумать над разбиением
void *pripole() {
	unsigned long long ex = 2500,	// количество разбиений системы по направлению X
		      ey = 2500;	// количество разбиений системы по направлению Y
	int i, j;
	double x = X / (2 * ex), dx = X / ex, y = Y / (2 * ey), dy = Y / ey;
	FILE *epole;

	epole = fopen("H.txt", "w");
	for (i = 0; i <= ex; i++)
	{
		for (j = 0; j <= ey; j++)
			fprintf(epole, "%.12f ", pole(x + dx * i,y + dy * j));
		fprintf(epole, "\n");
	}
	fflush(epole);
	fclose(epole);
}

// todo: подумать над разбиением
// todo: подумать учётом температуры
/* печатает файл с энергией взаимодействия вихрь-вихрь (E3) от расстояния */
void *pri3() {
	unsigned long long ex = 10000,	// количество разбиений системы по направлению X
		      ey = 2500;	// количество разбиений системы по направлению Y
	int i, j, k;
	double  x = X / (2 * ex), dx = X / ex, y = Y / (2 * ey), dy = Y / ey, r, xi, yi, xii, yii, v, s, f = 0, t;
	FILE *e3e3;

	e3e3 = fopen("E3.txt", "w");
	for (i = 0; i < ex; i++)
		fprintf(e3e3, "%.8f;%.16f;%.16f\n", (x + dx * i), e3(x + dx * i,lamb(1 / 11606)), L);

	fflush(e3e3);
	fclose(e3e3);
}


#define COPR	4
// использование:	pri(сколько будет печатей, что печатать, ...)
/* быстрый запуск необходимых функций для печати разного рода файлов с энергиями, разбиение для энергий берётся как от **ee */
void pri(int count, ...) {
	char ii[COPR];
	int i, c;
	va_list fa;
	pthread_t te[COPR];

	for (i = 0; i < COPR; i++)
		ii[i] = 0;

	va_start(fa, count);
	for (i = 0, c = va_arg(fa, int); i < count; i++, c = va_arg(fa, int))
		switch (c)
		{
			case '3':
				ii[0] = 1;
				pthread_create(&te[0], NULL, pri3, NULL);
				break;
			case 'p':
				ii[1] = 1;
				pthread_create(&te[1], NULL, primes, NULL);
				break;
			case 't':
				ii[2] = 1;
				pthread_create(&te[2], NULL, pritemp, NULL);
				break;
			case 'h':
				ii[3] = 1;
				pthread_create(&te[3], NULL, pripole, NULL);
				break;
			default:
				erro("такого ключа нет", 19); // erro 19
		}
	va_end(fa);

	for (i = 0; i < COPR; i++)
		if (ii[i]) pthread_join(te[i],NULL);
}




int main(int argc, char *argv[]) {
	char ss[512];
	int i, j, x, y, p, q;
	double f = 0, s, g, fs = 0;
	unsigned long long dk = 0;
	pthread_t th[3];


//	printf(">> %.16f\n>> %.16f\n", m1 * 400, e1(400,L0,0.00025,0.00050));
//	exit(0);

//	pri(1,'t');
//	exit(0);

	// генератор случайного
	srand(clock());


	// примеси
        pthread_create(&th[0], NULL, defect, NULL);

        pthread_join(th[0],NULL);


	// ввод из коммандной строки
	take(argc,argv);


	// печать всякой фигни
	pri(1,'p');


	// начальная энергия
	e = malloc(2 * sizeof(double));
	e[1] = e[0] = tote1(v[nm],nn[nm]) + tote2(v[nm],nn[nm]) + tote3(v[nm],nn[nm]) + tote4(v[nm],nn[nm]);

	me = malloc(2 * sizeof(double));
	me[1] = me[0] = e[0];


#if D < 0
	if (nn[nm] > 4) d = nn[nm];
	else d = DM;
#endif

	// основной цикл
S:	while (!kbhit())
	{
		pthread_create(&th[0], NULL, three, NULL);
		pthread_join(th[0],NULL);

		if (qqq) goto E;
	}
	switch (getch())
	{
		case 's':
			show();
			goto S;
		case 'q':
			// todo: сделать сохранение посчитоного до следующего програмного шага
			break;
		case 'p':
			// todo: сделать паузу, скушать твикс
			goto S;
		case 'i':
			// todo: сделать сохранение
			goto S;
		default:
			goto S;
	}

//E:	snprintf(ss, 512, "coordinat_%d.txt", rand() % 1000);
E:	snprintf(ss, 512, "coordinat.txt");
	save(ss);

	show();

	return 0;
}
