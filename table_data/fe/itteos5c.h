/*
 ITTEOS5.4, 15.10.2004 19:00

 (C) P.R.Levashov, K.V.Khishchenko, 2004
 The SOFTWARE PRODUCT and documentation are provided with RESTRICTED RIGHTS. 
 You can use or duplicate this product only by approval of the authors.
 If this software accidentally fell into your hands please contact
 the authors: "Pavel Levashov <pasha@ihed.ras.ru>, 
              "Konstantin Khishchenko <konst@ihed.ras.ru>
*/

/* 
   You can determine the values of parameters TABLES_MNS, TABLES_MCNT, TABLES_MNP
   by using the program params.for
*/
#define TABLES_MNS 1
#define TABLES_MCNT 4000
#define TABLES_MNP 257 


#define TABLES_NOV 12

extern void __stdcall eosini( double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ] );
extern void __stdcall  eostv( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  metatv( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  eospv( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  metapv( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  eosev( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  metaev( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  meltts( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  metameltts( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  melttl( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  metamelttl( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  bindtl( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  bindtg( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  sublts( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  subltg( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  spintl( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  spintg( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
extern void __stdcall  binodv( int* me, double tables[ TABLES_MNP ][ TABLES_MCNT ][ TABLES_MNS ], double eos[ TABLES_NOV ] );
