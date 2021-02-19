# 99_functions.R
# Feb 2021

process_results = function(in.results, max.category.num=3){

## construct probability of supporting vaccine donation from chains
names = names(in.results$sims.array[1,1,]) # names of all the stored results
all_results = NULL
for (countrynum in 1:n.country){
  index1 = grep(paste('^alpha\\[', countrynum, sep=''), names) # locations of the alpha's for this country
  indexi = grep(paste(',', 1 ,'\\]', sep=''), names) # intercept locations
  i = intersect(index1, indexi) # location of intercept for this country
  c1 = bugs.results$sims.array[,1,i] # chain 1  
  c2 = bugs.results$sims.array[,2,i] # chain 2
  frame = data.frame(countrynum = countrynum, num = 1, res=c(c1, c2), index=1:(length(c1)*2))
  all_results = bind_rows(all_results, frame)
  # other category levels
  for (num in 2:max.category.num){
    index2 = grep(paste(',', num ,'\\]', sep=''), names) # locations of these category levels
    u = intersect(index1, index2) # location of category (alpha + num)
    c1 = bugs.results$sims.array[,1,i] + bugs.results$sims.array[,1,u] # chain 1  
    c2 = bugs.results$sims.array[,2,i] + bugs.results$sims.array[,2,u] # chain 2
    frame = data.frame(countrynum = countrynum, num = num, res=c(c1, c2), index=1:(length(c1)*2))
    all_results = bind_rows(all_results, frame)
  }
}
# now add the overall means
indexi = grep('^mu.alpha\\[1\\]', names) # mean alpha for intercept
c1 = bugs.results$sims.array[,1,indexi] # chain 1  
c2 = bugs.results$sims.array[,2,indexi] # chain 2
frame = data.frame(countrynum = n.country+1, num = 1, res=c(c1, c2), index=1:(length(c1)*2)) # country number is + 1
all_results = bind_rows(all_results, frame)
for (num in 2:max.category.num){
  index2 = grep(paste('^mu.alpha\\[', num, '\\]', sep=''), names) # mean alpha for other questions
  c1 = bugs.results$sims.array[,1,indexi] + bugs.results$sims.array[,1,index2] # chain 1  
  c2 = bugs.results$sims.array[,2,indexi] + bugs.results$sims.array[,2,index2] # chain 2
  frame = data.frame(countrynum = n.country+1, num = num, res=c(c1, c2), index=1:(length(c1)*2))
  all_results = bind_rows(all_results, frame)
}
# 
all_results = mutate(all_results,
                     p = inv.logit(res), # inverse logit to get probability
                     overall = 1 + as.numeric(countrynum > n.country)) # binary for overall vs countries

return(all_results)

} # end of function
