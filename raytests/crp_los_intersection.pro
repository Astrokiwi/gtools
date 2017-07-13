function crp_los_intersection, part1, part2, part, n_match
; get particles that intersect the line-of-sight
;

  ;compile_opt idl2, hidden
  ;!EXCEPT = 0

  ; set defaults
  if n_elements(part) eq 0 then return,!NULL
  n_part = n_elements(part)
  idx = lindgen(n_part)
  jdx = where(part.idx eq part1.idx,m_match)
  if m_match gt 0 then idx = remove_value(jdx[0],idx)
  jdx = where(part.idx eq part2.idx,m_match)
  if m_match gt 0 then idx = remove_value(jdx[0],idx)
  n_part = n_elements(idx)
  if n_part eq 0 then return,!NULL

  ; setup line vector
  l_length = sqrt(total((part2.pos_cart_half-part1.pos_cart_half)^2.))
  l_norm = (part2.pos_cart_half-part1.pos_cart_half)/l_length
  l_ref = part1.pos_cart_half

  ; calculate intersection distances (in units of rcl)
  d_comp1 = rebin(l_ref,[3,n_part]) - part[idx].pos_cart_half
  d_comp2 = total(d_comp1*rebin(l_norm,[3,n_part]),1)
  if n_part gt 1 then begin
     d_comp2 = transpose(rebin(d_comp2,[n_part,3]))*rebin(l_norm,[3,n_part])
     d_los = sqrt(total((d_comp1-d_comp2)^2.,1))
     ; check if cloud is within interval of the two points
     t = -(-total(rebin(l_norm,[3,n_part])*part[idx].pos_cart_half,1) + replicate(total(l_norm*l_ref),n_part))
     ref_size = part1.size + abs(t/l_length)*(part2.size-part1.size)
  endif else begin
     d_comp2 = replicate(d_comp2,3) * l_norm
     d_los = sqrt(total((d_comp1-d_comp2)^2.))
     ; check if cloud is within interval of the two points
     t = -(-total(l_norm*part[idx].pos_cart_half) + total(l_norm*l_ref))
     ref_size = part1.size + abs(t/l_length)*(part2.size-part1.size)
     if ((t ge 0.) and (t le l_length) and (d_los lt ref_size+part[idx].size)) then return,part[idx] else return,!NULL
  endelse

  ; extract intersecting clouds and send back
  jdx = where((t ge 0.) and (t le l_length) and (d_los lt (ref_size+part[idx].size)),n_match)
  if n_match gt 0 then return, part[idx[jdx]] $
  else return, !NULL
  
end
