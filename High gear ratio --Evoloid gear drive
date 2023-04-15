enable_variable_cache = true;

-------------------------------------------------
--module/subject:Case study cyber physical production  system using AM 
--Topic :Involute_Gearing:Evoloid Drive with High Gear Ratio


--Creating a demonstartion model of High Ratio Evoloid Drive 
-------------------------------------------------





-------------------------------------------------
--Input parameters needed for gear design

z1=ui_number("Number of Teeth on pinion",4,1,6)
                      --input for number of teeth on pinion/driver,in an evoloid drive where number of teeths on the pinion less than 5 or 7 (The pinion/driver gear comprising less tooth contributes to the high possible transmission ratios)

z2=ui_number("Number of Teeth on external_gear",28,24,36)
                       --input for number of teeth on external_gear



m_t=ui_number("module",3,1,50)
                    -- Gear module

alpha_t=ui_numberBox("Pressure angle (Degrees)",20)
                   --pressure angle in degrees

Face_width=ui_number("Face Width(mm)",25,1,50)               --Gear Face Width

x_coef_2=ui_scalar("profile_shift_external_gear(mm)",-0.6,-2,2)
               --profile shift factor for external_gear(external_gear should have a negative profile shift coefficient),values as per table

x_coef_1=ui_scalar("profile_shift_Pinion(mm)",.822,-2,2)
               --profile shift factor for Pinion(pinion should have a positive profile shift coefficient),values as per table

h_ap_1=ui_scalar("addendummodifcoef_Pinion(mm)",0.659,-2,2)--addendum modification coefficient of pinion,values as per table

h_fp_1=ui_scalar("dedendummodifcoef_Pinion(mm)",1.1,-2,2)--dendendum modification coefficient of pinion,values as per table

h_ap_2=ui_scalar("addendummodifcoef_gear(mm)",1.1,-2,2)----addendum modification coefficient of gear,values as per table

h_fp_2=ui_scalar("dedendummodifcoef_gear(mm)",.66,-2,2)--dendendum modification coefficient of gear,values as per table

Rotation=math.floor(ui_scalar("evoloid drive rotation",0,0,360))
          --rotation in which the pinion gear rotating with the external gear)

helix_angle_value=ui_scalar("helix angle(Degrees)",45,20,45)
        -- helix angle of evoloid drive
-------------------------------------------------





-------------------------------------------------
--Function for generating the gear profiles

function gear(z,m_t,alpha_t_rad,x_coef,f_r,b,h_ap,h_fp)

local xy={}   --create the matrix => table, use local for local variables!

alpha_t_rad=alpha_t*math.pi/180    --pressure angle of gear in radian

d_ref=z*m_t  --reference diameter

r_ref=d_ref/2   --reference radius

d_b=d_ref*math.cos(alpha_t_rad)	--base diameter

r_b=d_b/2   --base radius
r_a=r_ref + m_t*(h_ap+x_coef)--Tip radius
r_f=r_ref - m_t*(h_fp-x_coef)--Root radius
f_r = r_b - r_f --------fillet radius	
h_a = m_t *h_ap--addendum
h_f = m_t *h_fp--dedendum	
--For any radius of the circle, below equation  can get the value of corresponding tooth thickness(modifying method)
S_ref=m_t*((math.pi/2)+ 2*x_coef* math.tan(alpha_t_rad))					-- tooth thickness of the reference circle

inv_alpha_t=math.tan(alpha_t_rad) - alpha_t_rad	  --where  alpha_t_rad is pressure angle

d_y = math.sqrt(math.pow(d_ref*math.sin(alpha_t_rad)- 2*(h_a-(m_t*x_coef)-h_f *(1-math.sin(alpha_t_rad))),2) + d_b * d_b)	
                --where d_y is the diameter of any circle

r_y=d_y/2  --radius of any circle

alpha_y=math.acos((d_ref*math.cos(alpha_t_rad))/d_y)

inv_alpha_y=math.tan(alpha_y) - alpha_y
  -- where alpha_y pressure angle in any circle

S_ty = d_y*((S_ref/d_ref)+inv_alpha_t - inv_alpha_y)+(x_coef*math.tan(alpha_y))

  --where S_ty is tooth thickness in any circle

 -- omega = (s_ty)/(0.5*d_y)  --angle swept along the form circle 

n_points=30   --for iterating the points

invo_angle= m_t*((math.pi/2) + 2*x_coef*math.tan(alpha_t_rad)) /r_ref +2*math.tan(alpha_t_rad) -2*alpha_t_rad
--To draw the involute between two circles w.r.t angle

function curve(r_1,r_2)
return math.sqrt(((r_2 * r_2) -(r_1 *r_1))/(r_1 *r_1))
end   
--function for drawing involute curve between two radii

curve_1 = curve(r_b,r_b)
curve_2 = curve(r_b,r_a)
--it is intended to start and stop the involute tooth between r_a and r_b.

--Fillet formation 
--Fillet formed at the root area of gear with slope of the line and with f_r(fillet radius) and r_b.(end point from the tooth profile above the base circle has the same common tangent with the start point from the tooth profile below the base circle)
--next step is to create parts under the base circle by finding the boundary of the same common tangent on the profile curve, that is, the position where the base circle curve intersects the profile curve.

--Firstly, get the common tangency’s slope.
function slope(contour) 
return ((contour[2].y - contour[1].y) / (contour[2].x - contour[1].x))	--slope=(y2-y1/x2-x1)	
end
--slope of line.

function involute(r_b,angle)			
return v(r_b *(math.sin(angle) - angle *math.cos(angle)), r_b *(math.cos(angle) + angle *math.sin(angle)))	
end	
--generate the involute at one end  w.r.t the 	radius of base circle and angle.this angle is defined in between two radii

function involute_mirror(involut_points)
return v(-involut_points.x,involut_points.y)
end  --mirroring the involute points generated above where y coordinates not changing

local points = {}  --create the matrix => table, use local for local variables!
for i = 1, n_points do
points[i] =involute(r_b,(curve_1 + (curve_2 - curve_1) * i / n_points))
end

m_s = slope(points)  --With the value slope, the next step is to obtain the  inverse slope.
inv_slope = math.atan(m_s)

line_parellel = {}	--create the matrix => table, use local for local variables

line_parellel[1] = v(points[1].x + f_r * math.cos(inv_slope + math.pi / 2),points[1].y + f_r * math.sin(inv_slope + math.pi / 2)) --A line parallel to the slope of line is formed so as to find the point to form the circle

d = (line_parellel[1].y - m_s * line_parellel[1].x) / math.sqrt(m_s * m_s + 1)
th1 = math.asin(d / (f_r + r_b)) + inv_slope
---- The point of intersection of the fillet circle and the base radius gives the fillet radius.

fillet_center = {} --create the matrix => table, use local for local variables
fillet_center[1] = v((f_r + r_f) * math.cos(th1), (f_r + r_f) * math.sin(th1))
th31 = 2 * math.pi + math.atan(fillet_center[1].y / fillet_center[1].x)
th33 = 3 * math.pi / 2 + inv_slope

function point_rotate(ang,points)
return v(math.cos(ang)*points.x + math.sin(ang)*points.y,math.cos(ang)*points.y - math.sin(ang)*points.x)
end --rotates the single tooth profile with respect to the origin

function circle(a,b,r,th)
return v(a+ r*math.cos(th),b+ r*math.sin(th))
end  --function for creating a circle 

--starting full gear profile including fillet
for i=1,z do 

for j=1,n_points do 
xy[#xy+1]=point_rotate(2*math.pi*i/z,circle(fillet_center[1].x,fillet_center[1].y,f_r,(th31 +(th33-th31) * j / n_points)))
end    --first fillet is formed 

curve_1 = curve(r_b,r_y)
curve_2 = curve(r_b,r_a)
--it is intended to start and stop the involute tooth between r_a and r_b

for j=1,n_points do 
xy[#xy+1]=point_rotate(2*math.pi*i/z,involute(r_b, (curve_1 +(curve_2-curve_1) *j / n_points)))
end      --first involute curve or the one side of involute points is formed ,it is done by using involute function to generate curve but inorder to generate the curve between two radii we had used point_rotate

--now  we have to mirror the involute points generated above where y coordinates not changing
for j=n_points,1,-1 do 
xy[#xy+1]=point_rotate(2*math.pi*i/z,point_rotate(invo_angle,involute_mirror(involute(r_b,(curve_1 +(curve_2-curve_1) *j / n_points)))))
end   
--the second involute curve has been obtained by the involute_mirror function  
for j=n_points,1,-1 do 
xy[#xy+1]=point_rotate(2*math.pi*i/z,point_rotate(invo_angle,involute_mirror(circle(fillet_center[1].x,fillet_center[1].y,f_r,(th31 +(th33-th31) * j / n_points)))))
end --second fillet is formed  using the involute_mirror function
end
xy[#xy+1]=xy[1]  --   -- close contour set end point to start point
return xy -- return value of function: table of vectors
end  
-------------------------------------------------






-------------------------------------------------
--the involute spur gear can be changed into helical gear using the extrude function below
function extrude(Contour, angle, dir_v, scale_v, z_steps)
   -- extrude a Contour to a shape by turning the contour to angle in z_steps
   -- extrude a Contour in a dircetion given by the vector dir_v
   -- extrude a Contour with a scaling factor given by vector scale_v 
   -- Contour: a table of vectors as a closed contour (start point and end point is the same)
   -- angle: roation angle of contour along z_steps in deg
   -- dir_v: vector(x,y,z) direction of extrusion
   -- sacle_v: vector(x,y,z) scaling factors for extrudion
   -- z_steps: number of steps for the shape, mostly z_steps=2 if angle is equal zero 
   local n_cont= #Contour
   local angle= angle/180*math.pi

   local Vertex= {}
   for j= 0,z_steps-1 do
      local phi= angle*j/(z_steps-1)
      local dir_vh= dir_v*j/(z_steps-1)
      local scale_vh= (scale_v - v(1,1,1))*(j/(z_steps-1)) + v(1,1,1)
      for i= 1,n_cont-1 do
          Vertex[i+j*n_cont]= v((dir_vh.x + scale_vh.x * (Contour[i].x*math.cos(phi) - Contour[i].y*math.sin(phi))),
                                (dir_vh.y + scale_vh.y * (Contour[i].x*math.sin(phi) + Contour[i].y*math.cos(phi))),
                                (dir_vh.z * scale_vh.z))
      end
      table.insert(Vertex,Vertex[1+j*n_cont])
   end

   local vertex_sum_1 = v(0,0,0)
   local vertex_sum_m = v(0,0,0)
   for i= 1,n_cont-1 do
      vertex_sum_1= vertex_sum_1 + Vertex[i]
      vertex_sum_m= vertex_sum_m + Vertex[i+n_cont*(z_steps-1)]
   end    
   table.insert(Vertex,vertex_sum_1/(n_cont-1)) --n_cont*m_cont + 1
   table.insert(Vertex,vertex_sum_m/(n_cont-1)) --n_cont*m_cont + 2


   Tri= {} -- !!! the index on table with Vertex starts with zero !!!
   local k= 1
   for j=0,z_steps-2 do
      for i= 0,n_cont-2 do
         Tri[k]=   v(i, i+1, i+n_cont) + v(1,1,1)*n_cont*j
         Tri[k+1]= v(i+1, i+n_cont+1, i+n_cont) + v(1,1,1)*n_cont*j
         k= k+2
      end
   end
   for i= 0,n_cont-2 do
      Tri[k]= v(i+1,i,n_cont*z_steps)
      k= k+1
   end
   for i= 0,n_cont-2 do
      Tri[k]= v(i+n_cont*(z_steps-1),i+1+n_cont*(z_steps-1),n_cont*z_steps+1)
      k= k+1
   end
   return(polyhedron(Vertex,Tri))
end
-------------------------------------------------




-------------------------------------------------
--Next the helix angle need to apply through a helix function into the extrude function under variable angle where (angle: roation angle of contour along z_steps )
function helicalangle(Face_width,z,m_t,helix_angle_value)  
d_p = z * m_t 
helix_angl = helix_angle_value*math.pi/180
return(math.tan(helix_angl)*Face_width*180/(d_p*math.pi))
end
-------------------------------------------------




-------------------------------------------------
--now just assign the main gear function and helix angle to the extrude function we get helical gear
Gear = gear(z2,m_t,alpha_t_rad,x_coef_2,f_r,Face_width,h_ap_2,h_fp_2)--assigning the gear tooth values to Gear
Full_gear = Gear --creating full gear tooth profile under Full_gear variable
helix_angle = helicalangle(Face_width,z2,m_t,helix_angle_value)  --assigning the helicalangle function to helix_angle
-------------------------------------------------





-------------------------------------------------
--to calculate the angle between the two sides of the tooth
function gearangle(m_t,x_coef,alpha_t,z)
alpha_t=alpha_t*math.pi/180
return(((math.pi*m_t/2) + 2*m_t*x_coef*math.tan(alpha_t))/(z*m_t/2) + 2*math.tan(alpha_t) - 2*alpha_t) 
end
-------------------------------------------------




-------------------------------------------------
--assigning the subfunctions to main function
gear_angle=gearangle(m_t,x_coef_2,alpha_t,z2)--assigning the gearangle function to gear_angle which will allow to calculate the angle between the two sides of the tooth
rotation_1=rotate(0,0,gear_angle*90/math.pi)
rotation_2=rotate(0,0,Rotation) --to rotate  gears in z axis with rotation_2
distance= m_t * 2
value_1=11.9
slot=translate(0,-value_1,0)*ccube(distance,distance*2,Face_width)--creating a slot for the shaft for the gear
function circle(r)
local x, y = 0, 0
local XY={}
for i = 1, 360 do
local angle = i * math.pi / 180
XY[i] = v(x + r * math.cos( angle ), y + r * math.sin( angle ))
end  
return XY
end--function for creating a circle 
Gear_cylinder=extrude(circle(distance),0,v(0,0,Face_width), v(1,1,1),20) --function for creating a cylinder
Gear_cylinder_diff=difference(Gear_cylinder,slot)--subtracting the slot(cube) from Gear_cylinder
Evoloid_gear = extrude(Full_gear ,helix_angle, v(0,0,-Face_width), v(1,1,1), 20) --extruding the helical gear 
Gear_cylinder_full=difference(Evoloid_gear,Gear_cylinder_diff)--subtracting the slot(cube) from Evoloid_gear
Evoloid_gear=rotation_2*difference(rotation_1*Gear_cylinder_full,translate(0,0,-Face_width)*Gear_cylinder_diff)
-------------------------------------------------




-------------------------------------------------
--text need to write according to this way
f = font(Path .. 'LiberationMono-Regular.ttf')
text = translate(15,20,0)*translate(0,0,0)*translate(0,-d_p/2,0)*translate(-15,-15,24)*scale(.15,.5,1) * f:str('Gear_Mechanism_Using_AM', 3)
text1 = translate(10,56,0)*translate(0,0,0)*translate(0,-d_p/2,0)*translate(-15,-15,24)*scale(.15,.5,1) * f:str('supervisor:Prof:Dr.stefan scherbath', 3)
text=union(text,text1)
-------------------------------------------------




-------------------------------------------------
--emit big gear
d=translate(10,26,0)*translate(0,0,25)*translate(0,-d_p/2,0)*Evoloid_gear
Evoloid_gear=difference(d,text)
emit(Evoloid_gear,128)
--creating the big  gear which will then mate with the pinion 
-------------------------------------------------






-------------------------------------------------
 emit(translate(10,-20,35)*rotate(90,0,0)*cylinder(3,15),128)
 emit(translate(10,3,35)*rotate(90,0,0)*cylinder(3,15),128)
-------------------------------------------------



-------------------------------------------------
--emitting the big gear shaft
shaft=cylinder(distance,Face_width+35)--creating a shaft
cube=translate(0,-value_1,0)*ccube(distance,distance*2,Face_width+5)--creating a cube
s1=translate(10,26,5)*translate(0,-d_p/2,-Face_width)*rotation_2*difference(shaft,cube)
s2=translate(10,-5,-11.8)*rotate(90,0,0)*cylinder(3.4,26)
s3=difference(s1,s2)
emit(s3,128)--emitting big gear shaft with hole downside of shaft for screw bolting
-------------------------------------------------




-------------------------------------------------
Gear = gear(z1,m_t,alpha_t_rad,x_coef_1,f_r,Face_width,h_ap_1,h_fp_1)--assigning the gear tooth values to Gear
Full_gear = Gear --creating full gear tooth profile under Full_gear variable
helix_angle = helicalangle(Face_width,z1,m_t,-helix_angle_value)  --assigning the helicalangle function to helix_angle
gear_angle=gearangle(m_t,x_coef_1,alpha_t,z1)--assigning the gearangle function to gear_angle which will allow to calculate the angle between the two sides of the tooth
rotation_1=rotate(0,0,-180 -(360/z1 - gear_angle*180/math.pi)/2)
rotation_2=rotate(0,0,-Rotation*z2/z1) --to rotate  gears in z axis with rotation_2
Evoloid_gear = extrude(Full_gear ,helix_angle, v(0,0,-Face_width), v(1,1,1), 20)--extruding the helical gear 
Evoloid_gear_pinion=rotation_2*(rotation_1*Evoloid_gear)
-------------------------------------------------


-------------------------------------------------
emit(translate(10,26,25)*translate(0,d_p/1.72,0)*Evoloid_gear_pinion,9)
--creating the first pinion gear which will mate with the external gear
-------------------------------------------------


-------------------------------------------------
emit(translate(10,26,25)*translate(0,-m_t*z2-(d_p/1.72),0)*Evoloid_gear_pinion,9)
--creating the second pinion gear which will mate with the external gear
-------------------------------------------------







-------------------------------------------------
--emitting the shaft for the pinion which has hole for M3 x 12mm bolt
s12=translate(10,26,58)*translate(0,-1.16*m_t*z2,0)*translate(0,d_p/2,-2*Face_width)*rotation_2*cylinder(.35*d_p,v(0,0,-30),v(0,0,-5))
s22=translate(10,-50,-11.8)*rotate(90,0,0)*cylinder(3.2,26)
s23=difference(s12,s22)
emit(s23,9)--emitting the shaft for the pinion which has hole for M3 x 12mm bolt
-------------------------------------------------




-------------------------------------------------
--the handle for small pinion 
a=cylinder(2.8,20)
b=cylinder(7,5)
c=cylinder(9,2)
d=translate(0,0,18)*difference(c,b)
e=union(a,d)
f=translate(0,0,19)*rotate(90,0,0)*cylinder(1,8)
g=union(e,f)
h=translate(0,8,19)*rotate(90,0,0)*cylinder(1,8)
i=union(h,g)
j=translate(0,0,19)*rotate(90,0,90)*cylinder(1,8)
k=union(j,i)
l=translate(-8,0,19)*rotate(90,0,90)*cylinder(1,8)
m1=union(l,k)
emit(translate(10,33,21)*m1,9)
--the handle calculations for the small pinion 
-------------------------------------------------





-------------------------------------------------
--emitting the shaft for the pinion which has hole for M3 x 12mm bolt
s121=translate(10,27,58)*translate(0,d_p/2,-2*Face_width)*rotation_2*cylinder(.35*d_p,v(0,0,-30),v(0,0,-5))
s22=translate(10,50,-11.8)*rotate(90,0,0)*cylinder(3.2,26)
s23=difference(s121,s22)
emit(s23,9)
-------------------------------------------------




-------------------------------------------------
print("Gear Ratio:" .. z2/z1)--the gear ratio for the evoloid drive
-------------------------------------------------




-------------------------------------------------
--creating base plate with the holes for M3 X 20mm bolt.
b=translate(8,-25,-32.9)*ccube(148,210,5)--base plate
c1=rotate(90,0,0)*rotate(90,0,0)*translate(60,-60,20)*ccylinder(3.4,51,20)--first hole for M3 X 20mm bolt to connect both plate and wooden plate
c2=rotate(90,0,0)*rotate(90,0,0)*translate(-40,-60,20)*ccylinder(3.4,51,20)--second hole for M3 X 20mm bolt to connect both plate and wooden plate
c3=rotate(90,0,0)*rotate(90,0,0)*translate(-40,110,20)*ccylinder(3.4,51,20)--third hole for M3 X 20mm bolt to connect both plate and wooden plate
c4=rotate(90,0,0)*rotate(90,0,0)*translate(60,110,20)*ccylinder(3.4,51,20)--fourth hole for M3 X 20mm bolt to connect both plate and wooden plate
c5=union(c1,c2)
c6=union(c5,c3)
c7=union(c6,c4)
emit(difference(b,c7))--creating base plate with the holes for M3 X 20mm bolt.
-------------------------------------------------





-------------------------------------------------
--creating first C-channel support
--first structure support for first pinion 
b1=translate(0,-100,0)*translate(0,80,-5)*translate(8,-48,2)*rotate(90,0,0)*ccube(80,5,30)--top surface first structure support for first pinion 
b2_1=translate(10,-65.4,0)*ccylinder((.35*d_p)+.4,51,20)
emit(difference(b1,b2_1))
b2=translate(0,-100,0)*translate(0,80,-48)*translate(45,-48,30)*rotate(0,90,90)*ccube(33,5,30)--sidewall for first structure support for first pinion 
emit(b2)
b3=translate(-74.5,-100,0)*translate(0,80,-48)*translate(45,-48,30)*rotate(0,90,90)*ccube(33,5,30)--sidewall for first structure support for first pinion 
emit(b3)
-------------------------------------------------






-------------------------------------------------
--creating second C-channel support
--second structure support for big gear
b11=translate(0,50,0)*translate(0,-100,0)*translate(0,80,-5)*translate(8,-48,2)*rotate(90,0,0)*ccube(80,5,30)--top surface second structure support for big gear
b2_1=translate(10,-16,0)*ccylinder(distance+.4,51,20)
emit(difference(b11,b2_1))
b21=translate(0,50,0)*translate(0,-100,0)*translate(0,80,-48)*translate(45,-48,30)*rotate(0,90,90)*ccube(33,5,30)--side wall second structure support for big gear
emit(b21)
b31=translate(0,50,0)*translate(-74.5,-100,0)*translate(0,80,-48)*translate(45,-48,30)*rotate(0,90,90)*ccube(33,5,30)--side wall second structure support for big gear
emit(b31)
-------------------------------------------------








-------------------------------------------------
--creating third C-channel support
--third structure support for second pinion
b41=translate(0,100,0)*translate(0,-100,0)*translate(0,80,-5)*translate(8,-48,2)*rotate(90,0,0)*ccube(80,5,30)--top surface first structure support for second pinion
b2_1=translate(10,33,0)*ccylinder((.35*d_p)+.4,51,20)
emit(difference(b41,b2_1))
b51=translate(0,100,0)*translate(0,-100,0)*translate(0,80,-48)*translate(45,-48,30)*rotate(0,90,90)*ccube(33,5,30)--sidewall for first structure support for second pinion
emit(b51)
b61=translate(0,100,0)*translate(-74.5,-100,0)*translate(0,80,-48)*translate(45,-48,30)*rotate(0,90,90)*ccube(33,5,30)--sidewall for first structure support for second pinion
emit(b61)
-------------------------------------------------
