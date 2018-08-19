//See Stam03 "Real-Time Fluid Dynamics for Games"

//vars and definations
var c=document.getElementById("mainCanvas");
var N=200,size=202*202;
var radius=200/N;
var radiusInt = parseInt(radius);
var dt=0.1, diff=0.0001, visc=0;
var force=3, source=200;
var mx,my;
var u=new Array(66*66), v=new Array(66*66), u_prev=new Array(66*66), v_prev=new Array(66*66);		
var dens=new Array(66*66), dens_prev=new Array(66*66);	
var f=0;						//frame account
var randomForceIsOn = false;
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

//clear data and canvas for first start, or reset them
function clearCanvas(){
	ctx.clearRect(0,0,200,200);
}
function clearData(){
	for (var i=0 ; i<size ; i++ ) {
		u[i] = 0;
		v[i] = 0;
		dens[i] = 0;
		dens_prev[i] = 0;
		u_prev[i] = 0;
		v_prev[i] = 0;
	}	
}
function reSet(){
	clearData();
	clearCanvas();
}

//the whole canvas was devided into N*N grid, a "dot" here means a 
//single cell size radius*radius, where radius = canvasWith / N
var ctx=c.getContext("2d");
function drawDot(x,y,color){
	if(color>255) color=255;
	if(color<0)   color=0;
	ctx.fillStyle="rgb("+color+","+color+","+color+")";
	ctx.fillRect(x,y,radiusInt+1,radiusInt+1);
}
//add a 10*10 source rect to the canvas
function add_source_rect(i,j){
	if(i<10||i>190||j<10||j>190) return;
	var x0=(i-9.5)*radius,x=(i+9.5)*radius;
	var y0=(j-9.5)*radius,y=(j+9.5)*radius;
	for(var i0=i-9;i0<i+9;i0++)
		for(var j0=j-9;j0<j+9;j0++)
			dens[IX(i0,j0)]+=source;
	for(;x0<x;x0++)
		for(;y0<y;y0++){
			drawDot(x,y,parseInt(dens[IX(i,j)]));
		}
}
//add source and force
c.onmousedown = function(e){
	mx = e.pageX;
	my = e.pageY;
}
c.onmouseup = function(e) {
	var old_mx = mx;
	var old_my = my;
	mx = e.pageX;
	my = e.pageY;
	if(mx<1||mx>199||my<1||my>199) return;
	if(old_mx<1||old_mx>199||old_my<1||old_my>199) return;
	if((mx == old_mx)&&(my == old_my)) add_source_rect(mx,my);
	u[IX(old_mx,old_my)] = force*(mx-old_mx);
	v[IX(old_mx,old_my)] = force*(my-old_my);
}
//draw density onto canvas according to dens[] array
function draw_density(dens,radius){
	var i,j,x,y;
	for(i=0;i<N;i++){
		x=(i-0.5)*radius;
		for(j=0;j<N;j++){
			y=(j-0.5)*radius;
			drawDot(x,y,parseInt(dens[IX(i,j)]))
		}
	}
}
//animate function
function animate(){
	requestAnimationFrame(animate);
	f++;
	dens_step ( dens, dens_prev, u, v, diff, dt );
	vel_step ( u, v, u_prev, v_prev, visc, dt );
	clearCanvas();
	draw_density(dens,radius);
	if(randomForceIsOn)	randomForce();
}
//add randomForce for mobile user who is unable to drag
function randomForce(){
	if(f%5==0){
		var rx = parseInt(Math.random()*200);
		var ry = parseInt(Math.random()*200);
		if(rx<1||rx>199||ry<1||ry>199) return;
		u[IX(rx,ry)] =100;
		v[IX(rx,ry)] =100;
		ctx.fillStyle = "red";
		ctx.fillRect(rx,ry,3,3);
	}
}
function randomForceFlag(){
	randomForceIsOn = !randomForceIsOn;
	randomForceStatus();
}
var rs= document.getElementById("randomForceFlag");
function randomForceStatus(){
	if(randomForceIsOn) rs.innerHTML="Random force is on";
	else rs.innerHTML="Random force is off";
}

//simulate start
function init(){
	ctx.fillStyle="white";
	ctx.fillRect(0,0,200,200);
	clearData();
	animate();
}
init();

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
	for ( k=0 ; k<25 ; k++ ) {
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
			x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];		
			if (x<0.5) x=0.5; if (x>N+0.5) x=N+0.5;	
			i0=parseInt(x); i1=i0+1;
			if (y<0.5) y=0.5; if (y>N+0.5) y=N+0.5;	
			j0=parseInt(y); j1=j0+1;										
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
{	
	//add_source ( x, x0, dt );						
	SWAP ( x0, x ); diffuse ( 0, x, x0, diff, dt );
	SWAP ( x0, x ); advect ( 0, x, x0, u, v, dt );
}
function vel_step ( u, v, u0, v0, visc, dt )
{
	add_source ( u, u0, dt ); add_source ( v, v0, dt );
	SWAP ( u0, u ); diffuse ( 1, u, u0, visc, dt );
	SWAP ( v0, v ); diffuse ( 2, v, v0, visc, dt );
	project ( u, v, u0, v0 );
	SWAP ( u0, u ); SWAP ( v0, v );
	advect ( 1, u, u0, u0, v0, dt ); advect ( 2, v, v0, u0, v0, dt );
	project ( u, v, u0, v0 );
}