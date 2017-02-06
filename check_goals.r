##jmd
##5.23.2011
##check_goals.r

check.goals <- function(sp, goals){ 
    ##goal prog
    #goal.rxns <- unique(bm.goals$rxn)
    fba.goal <- FBA(a=sp, fba.obj.rxns=NULL, goals=goals)
    v <- fba.goal$xopt
    goal.names <- union(intersect(goals, names(goals)), grep('GOAL-RXN', names(v), value=TRUE))
    vg <- v[goal.names]
    return(vg)
}
