/**
 * @file Matrix.js
 * @date 2013-06-26 13:37 PDT
 * @author Paul Reuter
 * @version 0.9.0
 *
 * @modifications
 * 0.9.0 - 2013-06-26 - Written but undergoing debugging.
 * 1.0.0 - 2013-06-26 - First release
 */



/**
 * Construct a new Matrix.
 * @see init
 */
function Matrix(/*Matrix|array|uint*/ m, /*mixed|uint*/ n, /*mixed*/ v) {
  this.rows = [];
  if( typeof(m) !== 'undefined' ) { 
    this.init(m,n,v);
  }
  return this;
};


Matrix.Vector = function(n,v) { 
  v = v || 0;
  var b = [];
  for(var i=0; i<n; i++) { 
    b[i] = v;
  }
  return b;
};


Matrix.NaN = Number.NaN;


Matrix.isNaN = function(/*number*/ value) {
  return (isNaN(value) || value===Matrix.NaN);
};


Matrix.isZero = function(/*number*/ value) {
  return (-0.0000001 < value && value < 0.0000001);
};


// TODO: What's the diff between map and each?
Matrix.prototype.map = function(/*function*/ cb) {
  var B = new Matrix(this);
  B.rows = B.rows.map(function(r) { return r.map(cb); });
  return B;
};


Matrix.prototype.each = function(/*function(v,i,j)*/ cb) {
  //var args = Array.prototype.slice.apply(arguments,1);  
  var A = new Matrix(this);
  var sA = A.size();
  for(var i=0,j=0,ni=sA[0],nj=sA[1]; i<ni; i++) {
    for(j=0; j<nj; j++) {
      A.rows[i][j] = cb.call(this,this.rows[i][j],i,j);
    }
  }
  return A;
};


Matrix.abs = function(/*Matrix|number*/ A) {
  return (Matrix.isMatrix(A)) ? A.map(Math.abs) : Math.abs(A);
};
Matrix.acos = function(/*Matrix|number*/ A) { 
  return (Matrix.isMatrix(A)) ? A.map(Math.acos) : Math.acos(A);
};
Matrix.asin = function(/*Matrix|number*/ A) { 
  return (Matrix.isMatrix(A)) ? A.map(Math.asin) : Math.asin(A);
};
Matrix.atan = function(/*Matrix|number*/ A) { 
  return (Matrix.isMatrix(A)) ? A.map(Math.atan) : Math.atan(A);
};
Matrix.atan2 = function(/*Matrix|number*/ A, /*Matrix|number*/ B) { 
  if( !Matrix.isMatrix(A) ) { 
    if( !Matrix.isMatrix(B) ) { 
      return new Matrix(1,1,Math.atan2(A,B));
    }
    A = new Matrix(B.size(),+A);
  } else if( !Matrix.isMatrix(B) ) { 
    B = new Matrix(A.size(),B);
  }
  return A.each(function(v,i,j) { return Math.atan2(v,B.rows[i][j]); });
};
Matrix.ceil = function(/*Matrix|number*/ A) { 
  return (Matrix.isMatrix(A)) ? A.map(Math.ceil) : Math.ceil(A);
};
Matrix.cos = function(/*Matrix|number*/ A) { 
  return (Matrix.isMatrix(A)) ? A.map(Math.cos) : Math.cos(A);
};
Matrix.exp = function(/*Matrix|number*/ A) { 
  return (Matrix.isMatrix(A)) ? A.map(Math.exp) : Math.exp(A);
};
Matrix.floor = function(/*Matrix|number*/ A) { 
  return (Matrix.isMatrix(A)) ? A.map(Math.floor) : Math.floor(A);
};
Matrix.log = function(/*Matrix|number*/ A) { 
  return (Matrix.isMatrix(A)) ? A.map(Math.log) : Math.log(A);
};
Matrix.max = function(/*Matrix|number*/ A) { 
  if( Matrix.isMatrix(A) ) {
    return Math.max.apply(Math,A.rows.slice(0).map(function(r) { return Math.max.apply(Math,r); }));
  }
  return Math.max.apply(Math,Array.prototype.slice.call(arguments));
};
Matrix.min = function(/*Matrix|number*/ A) { 
  if( Matrix.isMatrix(A) ) { 
    return Math.min.apply(Math,A.rows.slice(0).map(function(r) { return Math.min.apply(Math,r); }));
  }
  return Math.min.apply(Math,Array.prototype.slice.call(arguments));
};
Matrix.pow = function(/*Matrix|number*/ A, /*Matrix|number*/ B) {
  if( !Matrix.isMatrix(A) ) { 
    if( !Matrix.isMatrix(B) ) { 
      return new Matrix(1,1,Math.pow(A,B));
    }
    A = new Matrix(B.size(),+A);
  } else if( !Matrix.isMatrix(B) ) { 
    B = new Matrix(A.size(),B);
  }
  return A.each(function(v,i,j) { 
    return (Matrix.isZero(v)) ? 0 : Math.pow(v,B.rows[i][j]); 
  });
};
Matrix.random = function(/*Matrix|number*/ A) {
  return (Matrix.isMatrix(A)) ? A.map(Math.random) : Math.random(A);
};
Matrix.round = function(/*Matrix|number*/ A) { 
  return (Matrix.isMatrix(A)) ? A.map(Math.round) : Math.round(A);
};
Matrix.sin = function(/*Matrix|number*/ A) { 
  return (Matrix.isMatrix(A)) ? A.map(Math.sin) : Math.sin(A);
};
Matrix.sqrt = function(/*Matrix|number*/ A) { 
  return (Matrix.isMatrix(A)) ? A.map(Math.sqrt) : Math.sqrt(A);
};
Matrix.tan = function(/*Matrix|number*/ A) { 
  return (Matrix.isMatrix(A)) ? A.map(Math.tan) : Math.tan(A);
};


Matrix.dot = function(/*Matrix|number*/ A, /*Matrix|number*/ B) { 
  if( !Matrix.isMatrix(A) ) { 
    if( Matrix.isMatrix(B) ) { 
      return B.dot(A);
    }
    var C = new Matrix(1,1,+A);
    return C.dot(B);
  }
  return A.dot(B);
};

Matrix.div = function(/*Matrix|number*/ A, /*Matrix|number*/ B) { 
  if( !Matrix.isMatrix(A) ) { 
    if( Matrix.isMatrix(B) ) { 
      return B.div(A);
    }
    var C = new Matrix(1,1,+A);
    return C.div(B);
  }
  return A.div(B);
};

Matrix.add = function(/*Matrix|number*/ A, /*Matrix|number*/ B) {
  if( !Matrix.isMatrix(A) ) { 
    if( Matrix.isMatrix(B) ) { 
      return B.add(A);
    }
    var C = new Matrix(1,1,+A);
    return C.add(B);
  }
  return A.add(B);
};


Matrix.sub = function(/*Matrix|number*/ A, /*Matrix|number*/ B) {
  if( !Matrix.isMatrix(A) ) { 
    if( Matrix.isMatrix(B) ) { 
      return B.sub(A);
    }
    var C = new Matrix(1,1,+A);
    return C.sub(B);
  }
  return A.sub(B);
};


/**
 * Create and initialize a Matrix with zeros.
 * @param {uint|array} nr number of rows or result of Matrix.size() method.
 * @param {uint} nc number of cols if nr isn't an array.
 */
Matrix.zeros = function(/*uint|array*/ nr, /*uint*/ nc) { 
  if( typeof nc === 'undefined' && Array.isArray(nr) ) { 
    nc = nr[1];
    nr = nr[0];
  }
  return new Matrix(nr,nc,0);
};


/**
 * Create and initialize a Matrix with ones.
 * @param {uint|array} nr number of rows or result of Matrix.size() method.
 * @param {uint} nc number of cols if nr isn't an array.
 */
Matrix.ones = function(nr,nc) { 
  if( typeof nc === 'undefined' && Array.isArray(nr) ) { 
    nc = nr[1];
    nr = nr[0];
  }
  return new Matrix(nr,nc,1);
};


/**
 * Create and initialize a Matrix with zeros.
 * @param {uint} n number of rows and cols
 */
Matrix.eye = function(n) { 
  n = n || 1;
  var M = new Matrix(n,n,0);
  for(var i=0; i<n; i++) { 
    M.set(i,i,1);
  }
  return M;
};


/**
 * Test if object is a Matrix
 */
Matrix.isMatrix = function(/*Matrix*/ B) {
  return (B instanceof Matrix);
};


Matrix.prototype.replace = function(f,r) {
  var A = new Matrix(A);
  var sA = A.size();
  for(var i=0,j=0,ni=sA[0],nj=sA[1]; i<ni; i++) { 
    for(j=0; j<nj; j++) { 
      if( A.rows[i][j] === f ) { 
        A.rows[i][j] = r;
      }
    }
  }
  return A;
};


/**
 * Init an m-by-n matrix with value v.    m int, n int, v mixed.
 * OR: Init a size[] matrix with value n. m is array size 2, n is value.
 * OR: Init a matrix with another matrix. m is Matrix, n & v undef
 * OR: Init a matrix by a column vector.  m is array, n & v undef
 */
Matrix.prototype.init = function(/*Matrix|array|uint*/ m, /*uint|value*/ n, /*value*/ v) {

  if( typeof(v) === 'undefined' ) { 
    if( typeof(n) === 'undefined' ) { 
      if( Matrix.isMatrix(m) ) {
      // copy Matrix
        return this.copy(m);
      } else if( Array.isArray(m) ) { 
        if( m.length>0 && Array.isArray(m[0]) ) { 
        // convert array-of-array to Matrix
          return this.copy(m);
        } else { 
        // convert column-array to Matrix
          for(var i=0,n=m.length;i<n;i++) { 
            this.rows[i] = m[i];
          }
        }
      } else { 
        // if m and n aren't int: m should be a Matrix (to be copied). 
        return null;
      }
    } else {
      if( Array.isArray(m) && m.length == 2 ) {
      // m is size, n is value.
        v = v || n;
        n = m[1];
        m = m[0];
      } else if( parseInt(m) > 0 && parseInt(n) > 0 ) {
      // m is num-rows, n is num-cols, v is optional.
        v = v || 0;
      } else { 
        // m should be a size vector, n should be the fill value.
        return null;
      }
    }
  }
  if( m>0 && n>0 ) { 
    this.rows = [];
    for(var i=0,j=0; i<m; i++) {
      this.rows[i] = [];
      for(j=0; j<n; j++) { 
        this.rows[i][j] = v;
      }
    }
  }
  return this;
};


Matrix.prototype.copy = function(/*Matrix|array*/ other) {
  if( Matrix.isMatrix(other) ) { 
    other = other.get();
  }
  this.rows = other.slice(0).map(function(r) { return r.slice(0); });
  return this;
};


Matrix.prototype.add = function(/*Matrix|number*/ B) {
  var A = new Matrix(this);
  if( !Matrix.isMatrix(B) ) { 
    B = new Matrix(A.size(),+B);
  }

  var sA = this.size();
  var sB = B.size();
  if( sA[0]!=sB[0] || sA[1]!=sB[1] ) { 
    return null;
  }

  for(var i=0, ni=sA[0]; i<ni; i++) { 
    for(var j=0, nj=sA[1]; j<nj; j++) {
      A.rows[i][j] += B.rows[i][j];
    }
  }
  return A;
};


Matrix.prototype.sub = function(/*Matrix|number*/ B) {
  var A = new Matrix(this);
  if( !Matrix.isMatrix(B) ) { 
    B = new Matrix(A.size(),+B);
  }

  var sA = this.size();
  var sB = B.size();
  if( sA[0]!=sB[0] || sA[1]!=sB[1] ) { 
    return null;
  }

  for(var i=0, ni=sA[0]; i<ni; i++) { 
    for(var j=0, nj=sA[1]; j<nj; j++) {
      A.rows[i][j] -= B.rows[i][j];
    }
  }
  return A;
};


Matrix.prototype.find = function(/*optional bool function*/ test) {
  var sA = this.size(),
       B = Matrix.zeros(this.size());

  if( typeof(test) !== 'function' ) { 
    test = function(v) { return !!v; };
  }

  for(var r=0,c=0,nr=sA[0],nc=sA[1]; r<nr; r++) { 
    for(c=0; c<nc; c++) { 
      B.rows[r][c] = (test(this.rows[r][c])) ? 1 : 0;
    }
  }
  return B;
};



Matrix.prototype.T = function() { 
  var sA = this.size();
  var B = new Matrix(sA[1],sA[0],0);
  
  for(var i=0,ni=sA[0],nj=sA[1]; i<ni; i++) { 
    for(var j=0; j<nj; j++) { 
      B.rows[j][i] = this.rows[i][j];
    }
  }
  return B;
};


Matrix.prototype.cross = function(/*Matrix*/ B) { 
  var sA = this.size();
  var sB = B.size();
  if( sA[1] != sB[0] ) {
    console.log(sA + ' cannot cross '+sB);
    return null;
  }

  var C = Matrix.zeros(sA[0],sB[1]);
  for(var i=0, ni=sA[0]; i<ni; i++) { 
    for(var j=0, nj=sB[1]; j<nj; j++) { 
      for(var k=0, nk=sB[0]; k<nk; k++) { 
        C.rows[i][j] += this.rows[i][k] * B.rows[k][j];
      }
    }
  }
  return C;
};


Matrix.prototype.dot = function(/*Matrix*/ B) {
  if( !Matrix.isMatrix(B) ) { 
    B = new Matrix(this.size(),B);
  }

  var sA = this.size();
  var sB = B.size();
  if( sA[0]!=sB[0] || sA[1]!=sB[1] ) { 
    return null;
  }

  var C = new Matrix(this);
  for(var i=0, ni=sA[0]; i<ni; i++) { 
    for(var j=0, nj=sA[1]; j<nj; j++) { 
      C.rows[i][j] *= B.rows[i][j];
    }
  }
  return C;
};


Matrix.prototype.div = function(/*Matrix*/ B) {
  if( !Matrix.isMatrix(B) ) { 
    B = new Matrix(this.size(),B);
  }

  var sA = this.size();
  var sB = B.size();
  if( sA[0]!=sB[0] || sA[1]!=sB[1] ) { 
    return null;
  }

  var C = new Matrix(this);
  for(var i=0, ni=sA[0]; i<ni; i++) { 
    for(var j=0, nj=sA[1]; j<nj; j++) {
      C.rows[i][j] = (Matrix.isZero(B.rows[i][j])) ? Matrix.NaN : C.rows[i][j]/B.rows[i][j];
    }
  }
  return C;
};


Matrix.prototype.det = function() {
  var n = (this.size())[0],
    det = 0,
    sig = 1;

  var sA = this.size();
  if( sA[0]!=sA[1] ) { 
    return false;
  }
  if( n == 1 ) { 
    return this.rows[0][0];
  }
  if( n == 2 ) { 
    return this.rows[0][0]*this.rows[1][1] - this.rows[0][1]*this.rows[1][0];
  }

  for(var r in this.rows) {
    var col = this.rows[0][r];
    var B = new Matrix();
    for(var i=1; i<n; i++) { 
      var row = [];
      for(var j=0; j<n; j++) { 
        if( j!==col ) { 
          row.push(this.rows[i][j]);
        }
      }
      B.insertRow(row);
    }
    det += sig * col * B.det();
    sig *= -1;
  }
  return det;
};


Matrix.diag = function(v) {
  if( Matrix.isMatrix(v) ) { 
    return v.diag();
  }
  if( !Array.isArray(v) ) { 
    v = [v];
  }
  var n = v.length;
  var B = new Matrix(n,n,0);
  for(var i=0; i<n; i++) { 
    B.rows[i][i] = v[i];
  }
  return B;
};


Matrix.prototype.diag = function() { 
  var b = [];
  var sA = this.size();
  for(var i=0,ni=Math.min(sA[0],sA[1]); i<ni; i++) { 
    b.push(this.rows[i][i]);
  }
  return b;
};


Matrix.prototype.inv = function() {
  var sA = this.size();
  if( sA[0]!=sA[1] ) {
    return null;
  }
  var sA = this.size();
  var B = this.join(Matrix.eye(sA[0])).reduce();
  if( B !== 0 ) {
    return (B.split(sA[0]))[1];
  }
  return null;
};


Matrix.prototype.reduce = function() { 
  var B = new Matrix(this);
  var sB = this.size();
  
  for(var i=0,ni=sB[0],nj=sB[1]; i<ni; i++) { 
    var r = i, isZ = 1;
    
    while( r<ni ) { 
      if( !Matrix.isZero(B.rows[r][i]) ) { 
        isZ = 0;
        break;
      }
      r++;
    }
    if( isZ ) { 
      // inverse matrix doesn't exist
      return 0;
    }
    if( r > i ) { 
      var tmp = B.rows[i].slice(0);
      B.rows[i] = B.rows[r];
      B.rows[r] = tmp;
    }

    coeff = B.rows[i][i];
    for(var j=i; j<nj; j++) { 
      B.rows[i][j] /= coeff;
    }

    for(r=i+1; r<ni; r++) { 
      factor = B.rows[r][i];
      for(var j=0; j<nj; j++) { 
        B.rows[r][j] -= factor*B.rows[i][j];
      }
    }
  }
  
  for(var elim=ni-1; elim>0; elim--) { 
    for(var row=elim-1; row>=0; row--) { 
      var coeff = B.rows[row][elim];
      for(var j=0; j<nj; j++) { 
        B.rows[row][j] -= coeff*B.rows[elim][j];
      }
    }
  }

  return B;
};


Matrix.prototype.svd = function() {
// From Numerical Recipes
  var sA = this.size(), m = sA[0], n = sA[1],
      U = new Matrix(this),
      w = Matrix.Vector(n),
      V = new Matrix(n,n,0);

  var FMAX = function(a,b) { return (a>b) ? a : b; };
  var IMIN = function(a,b) { return (a<b) ? a : b; };
  var SIGN = function(a,b) { return (b<0.0) ? -Math.abs(a) : Math.abs(a); };
  var SQR = function(a) { return a*a; };
  var pythag = function(a,b) {
    var aa = Math.abs(a);
    var ab = Math.abs(b);
    if( aa>ab ) return aa*Math.sqrt(1+SQR(ab/aa));
    return (Matrix.isZero(aa)) ? 0 : ab*Math.sqrt(1+SQR(aa/ab));
  };

  var a = U.rows, // pointer to destination U (also input)
      v = V.rows; // pointer to destination V

  var flag,i,its,j,jj,k,l,mn;
  var anorm,c,f,g,h,s,scale,x,y,z;
  
  var rv1 = Matrix.Vector(1+n);
  g = scale = anorm = 0.0;
//  for(i=1; i<=n; i++) { 
  for(i=0; i<n; i++) { 
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;
//    if( i <= m ) { 
    if( i < m ) { 
//      for(k=i; k<=m; k++) { 
      for(k=i; k<m; k++) { 
        scale += Math.abs(a[k][i]);
      }
      if( scale ) { 
//        for(k=i; k<=m; k++) { 
        for(k=i; k<m; k++) { 
          a[k][i] /= scale;
          s += a[k][i]*a[k][i];
        }
        f = a[i][i];
        g = -SIGN(Math.sqrt(s),f);
        h = f*g - s;
        a[i][i] = f - g;
//        for(j=1; j<=n; j++) { 
        for(j=0; j<n; j++) { 
//          for(s=0.0,k=i; k<=m; k++) {
          for(s=0.0,k=i; k<m; k++) {
            s += a[k][i]*a[k][j];
          }
          f = s/h;
//          for(k=i; k<=m; k++) {
          for(k=i; k<m; k++) {
            a[k][j] += f*a[k][i];
          }
        }
//        for(k=i; k<=m; k++) {
        for(k=i; k<m; k++) {
          a[k][i] *= scale;
        }
      }
    }
    w[i] = scale * g;
    g = s = scale = 0.0;
//    if( i<=m && i!=n ) { 
    if( i<m && i!=n-1 ) { 
//      for(k=l; k<=n; k++) {
      for(k=l; k<n; k++) {
        scale += Math.abs(a[i][k]);
      }
      if( scale ) { 
//        for(k=l; k<=n; k++) { 
        for(k=l; k<n; k++) { 
          a[i][k] /= scale;
          s += a[i][k]*a[i][k];
        }
        f = a[i][l];
        g = -SIGN(Math.sqrt(s),f);
        h = f*g - s;
        a[i][l] = f-g;
//        for(k=l; k<=n; k++) {
        for(k=l; k<n; k++) {
          rv1[k] = a[i][k]/h;
        }
//        for(j=l; j<=m; j++) { 
        for(j=l; j<m; j++) { 
//          for(s=0.0,k=l; k<=n; k++) {
          for(s=0.0,k=l; k<n; k++) {
            s += a[j][k]*a[i][k];
          }
//          for(k=l; k<=n; k++) {
          for(k=l; k<n; k++) {
            a[j][k] += s*rv1[k];
          }
        }
//        for(k=l; k<=n; k++) {
        for(k=l; k<n; k++) {
          a[i][k] *= scale;
        }
      }
    }
    anorm = FMAX(anorm,(Math.abs(w[i])+Math.abs(rv1[i])));
  }

  // accumulation of right-hand transformations
//  for(i=n; i>=1; i--) { 
  for(i=n-1; i>=0; i--) { 
//    if( i<n ) { 
    if( i<n-1 ) { 
      if( g ) { 
//        for(j=l; j<=n; j++) { // double division to avoid underflow.
        for(j=l; j<n; j++) { // double division to avoid underflow.
          v[j][i] = (a[i][j]/a[i][l])/g;
        }
//        for(j=l; j<=n; j++) { 
        for(j=l; j<n; j++) { 
//          for(s=0.0,k=l; k<=n; k++) {
          for(s=0.0,k=l; k<n; k++) {
            s+= a[i][k]*v[k][j];
          }
//          for(k=l; k<=n; k++) {
          for(k=l; k<n; k++) {
            v[k][j] += s*v[k][i];
          }
        }
      }
//      for(j=l; j<=n; j++) { 
      for(j=l; j<n; j++) { 
        v[i][j] = v[j][i]=0.0;
      }
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }

  // Accumulation of left-hand transformations.
//  for(i=IMIN(m,n); i>=1; i--) {
  for(i=IMIN(m,n)-1; i>=0; i--) {
    l = i+1;
    g = w[i];
//    for(j=l; j<=n; j++) {
    for(j=l; j<n; j++) {
      a[i][j] = 0.0;
    }
    if( g ) { 
      g = 1.0/g;
//      for(j=l; j<=n; j++) { 
      for(j=l; j<n; j++) { 
//        for(s=0.0,k=l; k<=m; k++) {
        for(s=0.0,k=l; k<m; k++) {
          s+= a[k][i]*a[k][j];
        }
        f = (s/a[i][i])*g;
//        for(k=i; k<=m; k++) {
        for(k=i; k<m; k++) {
          a[k][j] += f*a[k][i];
        }
      }
//      for(j=i; j<=m; j++) {
      for(j=i; j<m; j++) {
        a[j][i] *= g;
      }
    } else { 
//      for(j=i; j<=m; j++) {
      for(j=i; j<m; j++) {
        a[j][i] = 0.0;
      }
    }
    ++a[i][i];
  }
  
//  for(k=n; k>=1; k--) { // Diagonalization of the bidiagonal form: Loop over
  for(k=n-1; k>=0; k--) { // Diagonalization of the bidiagonal form: Loop over
    for(its=1; its<=30; its++) {   // singular values, and over allowed iterations.
      flag = 1;
//      for(l=k; l>=1; l--) { // test for splitting.
      for(l=k; l>=0; l--) { // test for splitting.
        nm = l-1;   // Note that rv1[1] is always zero.
        if( (Math.abs(rv1[l])+anorm) == anorm ) { 
          flag = 0;
          break;
        }
        if( (Math.abs(w[nm])+anorm) == anorm ) break;
      }
      if( flag ) { 
        c = 0.0; // Cancellation of rv1[1], if l > 1.
        s = 1.0;
        for(i=l; i<=k; i++) { 
          f = s*rv1[i];
          rv1[i] = c*rv1[i];
          if( (Math.abs(f)+anorm) == anorm ) break;
          g = w[i];
          h = pythag(f,g);
          w[i] = h;
          h = 1.0/h;
          c = g*h;
          s = -f*h;
//          for(j=1; j<=m; j++) { 
          for(j=1; j<m; j++) { 
            y = a[j][nm];
            z = a[j][i];
            a[j][nm] = y*c + z*s;
            a[j][i] = z*c - y*s;
          }
        }
      }
      z = w[k];
      if( l==k ) { // Convergence.
        if( z<0.0 ) { // singular value is made nonnegative.
          w[k] = -z;
//          for(j=1; j<=n; j++) {
          for(j=0; j<n; j++) {
            v[j][k] = -v[j][k];
          }
        }
        break;
      }
      if(its==30) { console.log("no convergence in 30 sdv iterations"); }
      x = w[l];
      nm = k-1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g = pythag(f,1.0);
      f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c = s = 1.0; // Next QR transformation:
      for(j=l; j<=nm; j++) { 
        i = j+1;
        g = rv1[i];
        y = w[i];
        h = s*g;
        g = c*g;
        z = pythag(f,h);
        rv1[j] = z;
        c = f/z;
        s = h/z;
        f = x*c + g*s;
        g = g*c - x*s;
        h = y*s;
        y *= c;
//        for(jj=1; jj<=n; jj++) { 
        for(jj=0; jj<n; jj++) { 
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = x*c + z*s;
          v[jj][i] = z*c - x*s;
        }
        z = pythag(f,h);
        w[j] = z; // Rotation can be arbitrary if z=0.
        if( z ) { 
          z = 1.0/z;
          c = f*z;
          s = h*z;
        }
        f = c*g + s*y;
        x = c*y - s*g;
//        for(jj=1; jj<=m; jj++) { 
        for(jj=0; jj<m; jj++) { 
          y = a[jj][j];
          z = a[jj][i];
          a[jj][j] = y*c + z*s;
          a[jj][i] = z*c - y*s;
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }
  return [U,w,V];
};


Matrix.prototype.pinv = function() { 
  var uwv = this.svd();
  return uwv[2].cross(Matrix.pow(Matrix.diag(uwv[1]),-1).T()).cross(uwv[0].T());
};


Matrix.prototype.solve = function(b) { 
  // TODO: verify || optimize
  return (this.T().cross(this)).pinv().cross(this.T().cross(new Matrix([b])));
};


/**
 * Return all rows or specific cell.
 * @param uint ri (optional) Row index
 * @param uint ci (optional) Column index
 */ 
Matrix.prototype.get = function(/*uint*/ ri, /*uint*/ ci) {
  if( typeof ri === 'undefined' ) { 
    return this.rows;
  }
  var sA = this.size();
  if( typeof ci === 'undefined' ) {
    return (sA[0]>ri) ? this.rows[ri] : null;
  }
  return (sA[0]>ri && sA[1]>ci) ? this.rows[ri][ci] : null;
};


Matrix.prototype.set = function(/*uint*/ ri, /*uint*/ ci,/*value*/ v) {
  var sA = this.size();
  if( sA[0] <= ri || sA[1] <= ci ) { 
    return false;
  }
  this.rows[ri][ci] = v;
  return this;
};


Matrix.prototype.row = function(/*uint*/ rix) {
  return this.rows[rix].slice(0);
};


Matrix.prototype.col = function(/*uint*/ cix) {
  var arr = [];
  for(var i in this.rows) { 
    arr.push(this.rows[i][cix]);
  }
  return arr;
};


Matrix.prototype.toString = function() {
  return '('+this.size().join(',')+'): '+
    '['+this.rows.map(function(v) { return '['+v.join(',')+']'; })+']';
};


/**
 * @param {array} row An array to insert before ri
 * @param {uint} ri Row index to insert row before.
 */
Matrix.prototype.insertRow = function(/*array*/ row,/*uint*/ ri) {
  var sA = this.size();
  if( !Array.isArray(row) || Array.isArray(row[0]) ) { 
    return false;
  }
  if( typeof(ri) === 'undefined' || sA[0] < ri ) { 
    ri = sA[0];
  } else { 
    while( ri<0 ) { 
      ri += sA[0];
    }
  }
  this.rows.splice(ri,0,'TODO');
  this.rows[ri] = row;
  return this;
};


Matrix.prototype.removeRow = function(/*uint*/ ri) {
  var sA = this.size();
  if( typeof(ri) === 'undefined' || sA[0] < ri ) { 
    ri = sA[0] - 1 ;
  } else { 
    while( ri<0 ) { 
      ri += sA[0];
    }
  }
  this.rows.splice(ri,1);
  return this;
};


Matrix.prototype.join = function(/*Matrix*/ m) {
// TODO: support multiple Matricies as arguments.
  if( Matrix.isMatrix(m) ) { 
    m = m.get();
  }
  if( !Array.isArray(m) || m.length < 1 || !Array.isArray(m[0]) ) {
    return new Matrix(this);
  }
  var sA = this.size();
  var sB = [m.length,m[0].length];

  var C = new Matrix(Math.max(sA[0],sB[0]),sA[1]+sB[1],0);
  for(var i=0,ni=sA[0],nj=sA[1]; i<ni; i++) { 
    for(var j=0; j<nj; j++) { 
      C.rows[i][j] = this.rows[i][j];
    }
  }
  var nc = sA[1];
  for(var i=0,ni=sB[0],nj=sB[1]; i<ni; i++) { 
    for(var j=0; j<nj; j++) { 
      C.rows[i][nc+j] = m[i][j];
    }
  }
  return C;
};


Matrix.prototype.split = function(/*uint*/ cix) {
// TODO: support multiple indicies as arguments.
  var A = new Matrix(), B = new Matrix();
  var sA = this.size();
  for(var i=0,ni=sA[0]; i<ni; i++) {
    A.rows.push( this.rows[i].slice(0,cix) );
    B.rows.push( this.rows[i].slice(cix) );
  }
  return [A,B];
};


/**
 * @param {array} vect A column vector of values to insert before column ci.
 * @param {uint} ci Column index to insert values before.
 */
Matrix.prototype.insertColumn = function(/*arary*/ vect, /*uint*/ ci) {
  var sA = this.size();
  if( typeof(ci) === 'undefined' || sA[1] < ci ) { 
    ci = sA[1];
  } else { 
    while( ci<0 ) { 
      ci += sA[1];
    }
  }
  for(var i in this.rows) { 
    this.rows[i].splice(ci,vect[i]);
  }
  return this;
};


Matrix.prototype.removeColumn = function(ci) {
  var sA = this.size();
  if( typeof(ci) === 'undefined' || sA[1] < ci ) { 
    ci = sA[1]-1;
  } else { 
    while( ci<0 ) { 
      ci += sA[1];
    }
  }
  var v = [];
  for(var i in this.rows) { 
    v.push(this.rows[i].splice(ci,1));
  }
  return this;
};


/**
 * return array of coordinates of [row,col] of non-empty entries.
 */
Matrix.prototype.sparse = function() {
  var coords = [];
  var sA = this.size();

  for(var r=0,c=0,nr=sA[0],nc=sA[1]; r<nr; r++) { 
    for(c=0; c<nc; c++) {
      if( this.rows[r][c] ) { 
        coords.push([r,c]);
      }
    }
  }
  return coords;
}; 


/**
 * Return array of [num rows,num cols]
 */
Matrix.prototype.size = function() {     
  return [this.rows.length,(this.rows.length>0)?this.rows[0].length:0];
};









/** DEBUG **/
var A = new Matrix([ [1,0,0,0,2],[0,0,3,0,0],[0,0,0,0,0],[0,4,0,0,0] ]);
var B = Matrix.pow(A,3).add(42);
//console.log(A);
console.log(B);
process.exit();
var M = new Matrix([ [1,0,0,0,2],[0,0,3,0,0],[0,0,0,0,0],[0,4,0,0,0] ]);
var inv = M.pinv();
console.log(inv);
//var SVD = M.svd();
//console.log(SVD[0]);
//console.log(SVD[1]);
//console.log(SVD[2]);
//console.log(Matrix.diag(SVD[1]));
