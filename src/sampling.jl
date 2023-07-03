function getproposal(p::Prior, nparams)
  newparams = zeros(Float64, nparams)
  update_newparams!(newparams, p)
  return newparams
end

update_newparams!(newparams,p::Prior) = update_newparams!(newparams, 1, p.distribution...)

@inline function update_newparams!(newparams, i, x, y...)
  newparams[i] = rand(x)
  update_newparams!(newparams, i + 1, y...)
end
@inline function update_newparams!(newparams,i,x)
  newparams[i] = rand(x)
end

