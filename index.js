//var
var twoPI = 2*Math.PI;
var N=64,size=66*66;
var radius=200/N;
var dt=0.1, diff=0, visc=0;					//步长、扩散系数、粘度系数
var force=5, source=100;						//鼠标拖拽产生的力、amount of density that will be deposited
var dvel=true;								//当前绘制的是速度还是密度视图

var u=new Array(66*66), v=new Array(66*66), u_prev=new Array(66*66), v_prev=new Array(66*66);		
var dens=new Array(66*66), dens_prev=new Array(66*66);	
	
//defination
function IX(i,j){
	return i+(N+2)*j;
}
function SWAP(x0,x){
	var tmp=new Array(66*66);
	for(var i=0 ; i<size ; i++ ){
		tmp[i]=x0[i];
		x0[i]=x[i];
		x[i]=tmp[i];
	}
}
//init
var c=document.getElementById("mainCanvas");
var ctx=c.getContext("2d");
function drawDot(x,y,radius,r,a){
	ctx.beginPath();
	ctx.arc(x,y,radius,0,twoPI);
	ctx.fillStyle = "rgba("+r+","+r+","+r+","+ a +")";
	ctx.fill();
	ctx.closePath();
}

function initDot(){
	for (var i=0 ; i<size ; i++ ) {
		u[i] = 0;
		v[i] = 0;
		dens[i] = 0;
		dens_prev[i] = 0;
		u_prev[i] = 0;
		v_prev[i] = 0;
	}		
	console.log("!!!!",dens);
	u[IX(30,30)]=200;
	v[IX(30,30)]=200;
	for(var i=20;i<40;i++){
		x=(i-0.5)*radius;
		for(var j=20;j<40;j++){
			y=(j-0.5)*radius;
			drawDot(x,y,radius,0,1);
			dens[IX(i,j)]=100;
			console.log(i,j,IX(i,j),dens[IX(i,j)]);
		}
	}
	console.log("!!!!",dens);
}
initDot();
function animate(){
	requestAnimationFrame(animate);
	vel_step ( u, v, u_prev, v_prev, visc, dt );
	dens_step ( dens, dens_prev, u, v, diff, dt );
	draw_density(dens,radius);
}
function draw_density(dens,radius){
	var i,j,x,y;
	for(i=0;i<N;i++){
		x=(i-0.5)*radius;
		for(y=0;y<N;y++){
			y=(j-0.5)*radius;
			
			d00 = dens[IX(i,j)];
			d01 = dens[IX(i,j+1)];
			d10 = dens[IX(i+1,j)];
			d11 = dens[IX(i+1,j+1)];
			
			drawDot(x,y,radius,d00,1);
			drawDot(x,y+1,radius,d01,1);
			drawDot(x+1,y,radius,d10,1);
			drawDot(x+1,y+1,radius,d11,1);
		}
	}
	//console.log(dens);
}
// setInterval(function(){
	// animate();
// },50);
animate();
//solver
function add_source (x,s,dt)
{
	var i, size=(N+2)*(N+2);
	for ( i=0 ; i<size ; i++ ) x[i] += dt*s[i];
}
function set_bnd (b, x )
{
	var i;

	for ( i=1 ; i<=N ; i++ ) {
		x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
		x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
		x[IX(i,0  )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
		x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];
	}
	x[IX(0  ,0  )] = 0.5*(x[IX(1,0  )]+x[IX(0  ,1)]);
	x[IX(0  ,N+1)] = 0.5*(x[IX(1,N+1)]+x[IX(0  ,N)]);
	x[IX(N+1,0  )] = 0.5*(x[IX(N,0  )]+x[IX(N+1,1)]);
	x[IX(N+1,N+1)] = 0.5*(x[IX(N,N+1)]+x[IX(N+1,N)]);
}
function lin_solve ( b, x, x0, a, c )
{
	var i, j, k;
	for ( k=0 ; k<5 ; k++ ) {
		for(i=1;i<=N;i++){
			for(j=1;j<=N;j++){
				x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c;
			}
		}
		set_bnd ( b, x );
		
	}
}
function diffuse (b,x,x0,diff,dt)
{
	var a=dt*diff*N*N;
	lin_solve (b, x, x0, a, 1+4*a );
}
function advect (b,d,d0,u,v,dt)
{
	var i, j, i0, j0, i1, j1;
	var x, y, s0, t0, s1, t1, dt0;
	dt0 = dt*N;
	for(i=1;i<=N;i++)
		for(j=1;j<=N;j++){
			//线性回溯上一时刻的位置(x,y)
			x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];		
			//边界情况，并取到(x,y)的四个邻居
			if (x<0.5) x=0.5; if (x>N+0.5) x=N+0.5;	
			i0=parseInt(x); i1=i0+1;
			if (y<0.5) y=0.5; if (y>N+0.5) y=N+0.5;	
			j0=parseInt(y); j1=j0+1;										
			//四个邻居到点(x,y)线性插值，求得点(x,y)的值
			s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
			d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
						 s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
		}
	set_bnd ( b, d );
}
function project ( u,v,p,div )
{
	var i, j;

	for(i=1;i<=N;i++)
		for(j=1;j<=N;j++){
			div[IX(i,j)] = -0.5*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)])/N;
			p[IX(i,j)] = 0;
		}
	set_bnd ( 0, div ); set_bnd ( 0, p );

	lin_solve ( 0, p, div, 1, 4 );

	for(i=1;i<=N;i++)
		for(j=1;j<=N;j++){
			u[IX(i,j)] -= 0.5*N*(p[IX(i+1,j)]-p[IX(i-1,j)]);
			v[IX(i,j)] -= 0.5*N*(p[IX(i,j+1)]-p[IX(i,j-1)]);
		}
	set_bnd ( 1, u ); set_bnd ( 2, v );
}
function dens_step ( x, x0, u, v, diff, dt )
{	console.log(dens[1340],7);
	add_source ( x, x0, dt );
	console.log(dens[1340],7);
	SWAP ( x0, x ); 
	console.log(dens[0],8);
	diffuse ( 0, x, x0, diff, dt );
	console.log(dens[0],9);
	SWAP ( x0, x ); advect ( 0, x, x0, u, v, dt );
}

function vel_step ( u, v, u0, v0, visc, dt )
{
	add_source ( u, u0, dt ); add_source ( v, v0, dt );console.log(dens[0],1);
	SWAP ( u0, u ); diffuse ( 1, u, u0, visc, dt );console.log(dens[0],2);
	SWAP ( v0, v ); diffuse ( 2, v, v0, visc, dt );console.log(dens[0],3);
	project ( u, v, u0, v0 );console.log(dens[0],4);
	SWAP ( u0, u ); SWAP ( v0, v );
	advect ( 1, u, u0, u0, v0, dt ); advect ( 2, v, v0, u0, v0, dt );console.log(dens[0],5);
	project ( u, v, u0, v0 );console.log(dens[0],6);
}