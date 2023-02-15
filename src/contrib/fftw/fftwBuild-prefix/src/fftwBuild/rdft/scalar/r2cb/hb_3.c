/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Sat Jul 30 16:49:40 EDT 2016 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2hc.native -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -n 3 -dif -name hb_3 -include hb.h */

/*
 * This function contains 16 FP additions, 14 FP multiplications,
 * (or, 6 additions, 4 multiplications, 10 fused multiply/add),
 * 27 stack variables, 2 constants, and 12 memory accesses
 */
#include "hb.h"

static void hb_3(R *cr, R *ci, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 4); m < me; m = m + 1, cr = cr + ms, ci = ci - ms, W = W + 4, MAKE_VOLATILE_STRIDE(6, rs)) {
	       E Tk, Tj, Tn, Tl, Tm, To;
	       {
		    E T1, Td, T7, T8, T4, Tg, T2, T3;
		    T1 = cr[0];
		    T2 = cr[WS(rs, 1)];
		    T3 = ci[0];
		    Td = ci[WS(rs, 2)];
		    T7 = ci[WS(rs, 1)];
		    T8 = cr[WS(rs, 2)];
		    T4 = T2 + T3;
		    Tg = T2 - T3;
		    {
			 E T5, Tc, Tf, Ta, T9, Te, T6, Th, Ti, Tb;
			 T5 = W[0];
			 T9 = T7 + T8;
			 Te = T7 - T8;
			 cr[0] = T1 + T4;
			 T6 = FNMS(KP500000000, T4, T1);
			 Tc = W[1];
			 ci[0] = Td + Te;
			 Tf = FNMS(KP500000000, Te, Td);
			 Tk = FMA(KP866025403, T9, T6);
			 Ta = FNMS(KP866025403, T9, T6);
			 Tj = W[2];
			 Tn = FNMS(KP866025403, Tg, Tf);
			 Th = FMA(KP866025403, Tg, Tf);
			 Ti = Tc * Ta;
			 Tb = T5 * Ta;
			 Tl = Tj * Tk;
			 Tm = W[3];
			 ci[WS(rs, 1)] = FMA(T5, Th, Ti);
			 cr[WS(rs, 1)] = FNMS(Tc, Th, Tb);
		    }
	       }
	       cr[WS(rs, 2)] = FNMS(Tm, Tn, Tl);
	       To = Tm * Tk;
	       ci[WS(rs, 2)] = FMA(Tj, Tn, To);
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 3},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 3, "hb_3", twinstr, &GENUS, {6, 4, 10, 0} };

void X(codelet_hb_3) (planner *p) {
     X(khc2hc_register) (p, hb_3, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2hc.native -compact -variables 4 -pipeline-latency 4 -sign 1 -n 3 -dif -name hb_3 -include hb.h */

/*
 * This function contains 16 FP additions, 12 FP multiplications,
 * (or, 10 additions, 6 multiplications, 6 fused multiply/add),
 * 15 stack variables, 2 constants, and 12 memory accesses
 */
#include "hb.h"

static void hb_3(R *cr, R *ci, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * 4); m < me; m = m + 1, cr = cr + ms, ci = ci - ms, W = W + 4, MAKE_VOLATILE_STRIDE(6, rs)) {
	       E T1, T4, Ta, Te, T5, T8, Tb, Tf;
	       {
		    E T2, T3, T6, T7;
		    T1 = cr[0];
		    T2 = cr[WS(rs, 1)];
		    T3 = ci[0];
		    T4 = T2 + T3;
		    Ta = FNMS(KP500000000, T4, T1);
		    Te = KP866025403 * (T2 - T3);
		    T5 = ci[WS(rs, 2)];
		    T6 = ci[WS(rs, 1)];
		    T7 = cr[WS(rs, 2)];
		    T8 = T6 - T7;
		    Tb = KP866025403 * (T6 + T7);
		    Tf = FNMS(KP500000000, T8, T5);
	       }
	       cr[0] = T1 + T4;
	       ci[0] = T5 + T8;
	       {
		    E Tc, Tg, T9, Td;
		    Tc = Ta - Tb;
		    Tg = Te + Tf;
		    T9 = W[0];
		    Td = W[1];
		    cr[WS(rs, 1)] = FNMS(Td, Tg, T9 * Tc);
		    ci[WS(rs, 1)] = FMA(T9, Tg, Td * Tc);
	       }
	       {
		    E Ti, Tk, Th, Tj;
		    Ti = Ta + Tb;
		    Tk = Tf - Te;
		    Th = W[2];
		    Tj = W[3];
		    cr[WS(rs, 2)] = FNMS(Tj, Tk, Th * Ti);
		    ci[WS(rs, 2)] = FMA(Th, Tk, Tj * Ti);
	       }
	  }
     }
}

static const tw_instr twinstr[] = {
     {TW_FULL, 1, 3},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 3, "hb_3", twinstr, &GENUS, {10, 6, 6, 0} };

void X(codelet_hb_3) (planner *p) {
     X(khc2hc_register) (p, hb_3, &desc);
}
#endif				/* HAVE_FMA */
