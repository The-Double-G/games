var REST_DENS = 400
var GAS_CONST = 2000
var H = 16
var HSQ = H*H
var MASS = 5
var VISC = 100
var POLY6 = Math.sqrt(H)/(Math.PI*Math.pow(H,8))
var SPIKY_GRAD = -(Math.sqrt(H)+6) / (Math.PI*H**5)
var VISC_LAP = Math.sqrt(H)*10 / (Math.PI*H**5)
var EPS = H
var BOUND_DAMP = -0.5
var G
var DT = 0.0007

class Particle{
  constructor(x,y){
    this.pos = createVector(x,y)
    this.v = createVector()
    this.rho = 0
    this.f = createVector()
    this.p = 0
    this.radius = H/2
    this.scu = createVector()
  }
  show(){
    circle(this.pos.x,this.pos.y,H/2)
  }
}

function densityPressure(){
  for (var p of particles){
    p.rho = 0
    for (var pj of particles){
      var rij = p5.Vector.sub(pj.pos,p.pos)
      var r2 = rij.magSq()
      if (r2 < HSQ){
        p.rho += MASS * POLY6 * pow(HSQ-r2,3)
      }
    }
    p.p = GAS_CONST*(p.rho-REST_DENS)
  }
}

function forces(){
  for (var p of particles){
    var fpress = createVector()
    var fvisc = createVector()
    for (var pj of particles){
      if (pj == p){
        continue
      }
      var rij = p5.Vector.sub(pj.pos,p.pos)
      var r = rij.mag()
      
      if (r < H){
        var a = p5.Vector.mult(rij.normalize(),-MASS*(pj.p+p.p)*SPIKY_GRAD*pow(H-r,3)/(2*pj.rho))
        
        
        var con = VISC*MASS*VISC_LAP*(H-r)/(pj.rho)
        fpress.add(a)
        fvisc.add(p5.Vector.mult(p5.Vector.sub(pj.v,p.v),con  ))
      }
    }
    var fgrav = p5.Vector.mult(G,MASS/p.rho)
    var fo = createVector()
    fo.add(fpress)
    fo.add(fvisc)
    fo.add(fgrav)
    p.f=fo
  }
}

function integrate(){
  for (var p of particles){
    p.v.add(p5.Vector.mult(p.f,DT/p.rho))
    p.pos.add(p5.Vector.mult(p.v,DT))
    if (p.pos.y+H>height){
      p.v.y*=BOUND_DAMP
      p.pos.y=height-EPS
    }
    if (p.pos.x+EPS>width){
      p.v.x*=BOUND_DAMP
      p.pos.x=width-EPS
    }
    if (p.pos.x<EPS){
      p.v.x*=BOUND_DAMP
      p.pos.x=EPS
    }
    if (!meta){
      p.show()
    }
  }
}