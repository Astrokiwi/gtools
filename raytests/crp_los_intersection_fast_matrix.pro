function crp_los_intersection_fast_matrix, ref_part, part, params
; get particles that intersect the line-of-sight
;
  
  ;compile_opt idl2, hidden
  ;!EXCEPT = 0

  ; set defaults
  if n_elements(part) eq 0 then return,!NULL
  n_part = n_elements(part)
  idx = where(part.idx ne ref_part.idx,n_idx)
  stoch = params.comp.stochastical
  idm = dblarr(n_idx) + 1.
  idm3 = dblarr(3) + 1

  res = dblarr(n_part,7)

  tf_mat = part[idx].pos_cart_half - (ref_part.pos_cart_half # idm)
  tf_phi = transpose(atan(tf_mat[1,*],tf_mat[0,*]))
  tf_theta = transpose(atan(tf_mat[2,*],sqrt(tf_mat[0,*]^2+tf_mat[1,*]^2)))
  tf_dis = sqrt(total(tf_mat^2,1))
  
  res[idx,2] = tf_dis
  rho = (0.5 - tf_dis/(2.*sqrt(tf_dis^2+part[idx].size^2)))
  res[idx,6] = rho
  rho = sqrt(temporary(4.*!PI*rho)/!PI)
  res[idx,3:5] = transpose(tf_mat[0:2,*])/(res[idx,2] # idm3)

  ang_dif = (sin(tf_theta) # sin(tf_theta)) + (cos(tf_theta) # cos(tf_theta)) * cos((tf_phi # idm) - (idm # tf_phi))
  ang_dif = acos(temporary(ang_dif))

  rhoc = (rho # idm)
  for i=0,n_idx-1 do rhoc[*,i] += rho[i]
  a_jdx = ((ang_dif lt rhoc) and (tf_dis # idm) lt (idm # tf_dis))
  for i=0,n_idx-1 do a_jdx[i,i] = 0
  
  res[idx,0] = total(a_jdx,1)
  
  kdx = where(res[idx,0] gt 0,k_match)
  if k_match gt 0 then begin
     for i=0,k_match-1 do begin
        ldx = where(a_jdx[*,kdx[i]] eq 1)
        if stoch eq 1 then begin
           f_eff = total(crp_poisson(res[idx[kdx[i]],0],/norm)*exp(-mean(part[idx[ldx]].tau)*findgen(res[idx[kdx[i]],0]+1)))
           res[idx[kdx[i]],1] = -ALOG(f_eff)
        endif else res[idx[kdx[i]],1] = total(part[idx[ldx]].tau)/sqrt(2.)
     endfor
  endif

  return, res
end
