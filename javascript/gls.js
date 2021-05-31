/* gls.js
The generalised Lomb-Scargle periodogram
                     (Zechmeister & Kuerster, 2009, A&A, 496, 577)

Author: M. Zechmeister (Institut for Astrophysics Goettingen)
Version: 2019-03-28 (1.00)

gls = GLS(t_data, y_data, kwargs)
kwargs = {e_y:e_y, fbeg:fbeg, fend:fend, ofac:ofac, ls:ls}
*/

function dot(x, y) {
   var i=x.length, sum=0.;
   while(i--) sum += x[i] * y[i];
   return sum;
};

function add(x, a) {
   var i=x.length, xa=[];
   while(i--) xa[i] = x[i] + a;
   return xa;
};

function mul(x, y) {
   var i=y.length, xy=[];
   if (x.length==i) {
      while(i--) xy[i] = x[i] * y[i];
   } else {
      while(i--) xy[i] = x * y[i];
   };
   return xy;
};


function GLS(t_data, y_data, kwargs) {
   var k, i, nt, t, tmin, tbase, y, w=[], wy, wsum=0.;
   var df, delta_t, fbeg, fend, nf;
   var f, p, kbest=0;
   var kwargs = kwargs || {};
   var ofac = kwargs.ofac || 20.;
   var twopi = 2.0 * Math.PI;

   tmin = Math.min.apply(null, t_data);
   t = add(t_data, -tmin);
   tbase = Math.max.apply(null, t);

   i = nt = t.length;
   while(i--) wsum += w[i]= (kwargs.e_y? 1./kwargs.e_y[i]/kwargs.e_y[i] : 1.);

   // normalize weights, now "wsum=1"
   i = nt;
   while(i--) w[i] /= wsum;

   ymean = dot(w, y_data);
   y = add(y_data, -ymean);
   wy = mul(w, y);
   YY = dot(wy, y);    // Eq.(10), variance for the zero mean data

   df = 1. / tbase / ofac;       // frequency sampling depends on the time span, default for start frequency
   delta_t = tbase / (nt-1);     // mean time sampling
   fbeg = parseFloat(kwargs.fbeg || df);
   fend = parseFloat(kwargs.fend || 0.5/delta_t);
   nf = Math.floor((fend-fbeg)/df)+1;  // size of frequency grid

   f = new Array(nf);
   p = new Array(nf);

   for(k=0; k<nf; k++) {
      f[k] = fbeg + k * df;
      omega = twopi * f[k];
      C = 0.;
      S = 0.;
      YC = 0.;
      YS = 0.;
      CC = 0.;
      SS = 0.;
      CS = 0.;
      for(i=0; i<nt; i++) {
         wi = w[i];
         cosx = Math.cos(omega * t[i]);
         sinx = Math.sin(omega * t[i]);
         if (!kwargs.ls) {
         C += wi * cosx;         // Eq.(8)
         S += wi * sinx;         // Eq.(9)
         };
         YC += wy[i] * cosx;     // Eq.(11) simplifies, since Y should be 0 (mean was subtracted)
         YS += wy[i] * sinx;     // Eq.(12)   -"-
         CC += wi * cosx * cosx; // Eq.(13)
         CS += wi * cosx * sinx; // Eq.(15)
      };
      SS = 1. - CC;           // Eq.(14)
      SS -= S * S;            // Eq.(14)
      CC -= C * C;            // Eq.(13)
      CS -= C * S;            // Eq.(15)
      D = CC*SS - CS*CS;      // Eq.(6)

      A = (YC*SS-YS*CS) / D
      B = (YS*CC-YC*CS) / D
      off = -A*C-B*S

      p[k] = (SS*YC*YC/D + CC*YS*YS/D - 2.*CS*YC*YS/D) / YY;  // Eq.(5)
      if (p[k]>p[kbest]) {kbest = k; Abest = A; Bbest = B; Cbest = off+ymean}
   };
   return {p:p, f:f, k:kbest, Ak:Abest, Bk:Bbest, Ck:Cbest}
};
