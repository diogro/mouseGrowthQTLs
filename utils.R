makeSimData = function(data, marker, percent){
  a = findA(data, marker, percent)
  data + marker * a
}
markerVar = function(a, data, marker){
    data = mean(data) + residuals(lm(data ~ marker))
    za = data + marker * a
    n = length(data)
    V_pa = sum((za - mean(za))^2)/(n - 1)
    gen_means = tapply(za, marker, mean)
    gen_n = table(marker) / (n - 3)
    var_marker = gen_n %*% (gen_means - mean(za))^2
    (var_marker / V_pa)[1]
}
findA = function(data, marker, percent){
  V_p = var(data)
  a = sqrt((2*percent*V_p)/(1 - percent))
  return(a)
}

containsZero = function(x) ifelse(0 > x[1] & 0 < x[2], TRUE, FALSE)
find_CI = function(x, prob = 0.95){
    n = length(x)
    xs = sort(x)
    nint = floor(prob*n)
    lowest_int = abs(xs[n] - xs[1])
    for(i in 1:(n-nint)){
        current_int = abs(xs[i] - xs[i+nint])
        if(current_int <= lowest_int){
            lowest_int = current_int
            pos = i
        }
    }
    return(c(xs[pos], xs[pos+nint]))
}
