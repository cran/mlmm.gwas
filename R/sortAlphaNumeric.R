">.alphnum" = function(a,b){
    a = as.character(a)
    b = as.character(b)

    if(a==b) return(TRUE)

    x = c(a,b)
    for(i in 1:2){
        x=gsub("([0-9])([^0-9~])|([^0-9~])([0-9])","\\1\\3~\\2\\4", x)
    }
    x = strsplit(x, "~")
    a=x[[1]]
    b=x[[2]]

    for(i in 1:(max(length(a),length(b)))){
        if(a[i] != b[i]){
            if(is.na(b[i])) return(TRUE)
            if(is.na(a[i])) return(FALSE)
            if(suppressWarnings(!is.na(as.numeric(a[i])) & !is.na(as.numeric(b[i])))){
                ai = as.numeric(a[i])
                bi = as.numeric(b[i])
            }else{
                ai = a[i]
                bi = b[i]
            }
            return(ai > bi)
        }
    }

    return(TRUE)
}
"==.alphnum" = function(a,b){ as.character(a) == as.character(b) | !( a < b | b < a ) }
"[.alphnum" = function(x, i){
    class(x) = 'character'
    x = x[i]
    class(x) = "alphnum"
    x
}
sortAlphaNumeric = function(x){
    class(x) = "alphnum"
    x = sort(x)
    class(x) = "character"
    x
}
