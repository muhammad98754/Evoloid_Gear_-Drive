emit(translate(8,-25,-32.9)*ccube(148,210,5),10 )--base plate
--emit(translate(8,-128,2)*rotate(90,0,0)*ccube(148,74,5),10 )--sidewall

--emit(translate(0,206,0)*translate(8,-128,2)*rotate(90,0,0)*ccube(148,74,5),21 )--sidewall

--emit(translate(-64,-25,2)*rotate(0,90,0)*ccube(74,210,5),3 )--sidewall

--emit(translate(80,-25,2)*rotate(0,90,0)*ccube(74,210,5),1 )--sidewall

emit(translate(0,80,12)*translate(8,-128,2)*rotate(90,0,0)*ccube(80,10,30),21 )--sidewall

emit(translate(0,80,-10)*translate(40.5,-128,2)*rotate(90,0,0)*ccube(15,54,30),5 )--sidewall

emit(translate(0,80,-10)*translate(-38,-128,2)*rotate(90,0,0)*ccube(15,54,30),5 )--sidewall

--second one support

emit(translate(0,80,12)*translate(8,-48,2)*rotate(90,0,0)*ccube(80,10,30),21 )--sidewall

emit(translate(0,80,-10)*translate(40.5,-48,2)*rotate(90,0,0)*ccube(15,54,30),5 )--sidewall

emit(translate(0,80,-10)*translate(-38,-48,2)*rotate(90,0,0)*ccube(15,54,30),5 )--sidewall

function circle_offset(r,o)
  local tbl = {}
  local n = 128
  for i=1,n do 
   local a = 360*i/n
   tbl[i] = r * v(cos(a),sin(a),0) + o
  end
  return tbl
end

function bolt(head_sz,head_h)
  local head = intersection{
    cube(head_sz,30,head_h),
    rotate(0,0,60)*cube(head_sz,30,head_h),
    rotate(0,0,120)*cube(head_sz,30,head_h),
    }
  return head
end



function screw(ra,rs,len)
  local all_tbl={}
  local n = len*10
  local a
  for h = 1,n do
    a = h * 10
    all_tbl[h] = circle_offset(ra,v(rs*cos(a),rs*sin(a),h/10))
  end
  local head_sz = ra*2.5
  local head_h  = ra/1.2
  local head    = bolt(head_sz,head_h)
  local s = union{
    translate(0,0,head_h)*sections_extrude(all_tbl),
    translate(0,0,head_h)*translate(rs*cos(a),rs*sin(a),n/10)*cone(ra,ra/2,ra/2),
    head
    }
  return s
end

emit(rotate(90,0,90)*translate(-48,-5,-60)*screw(7,0.5,120))

b = translate(0,0,0)*difference( 
   bolt(20,10),
   translate(0,0,-10)*screw(7.3,0.5,20)
   )
emit(translate(55,-48,-5)*rotate(90,0,90)*b,21)

emit(rotate(90,0,90)*translate(33,-5,-60)*screw(7,0.5,120))

emit(translate(55,33,-5)*rotate(90,0,90)*b,21)
