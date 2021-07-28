function solution = newtonRaphson(prevSolution, prevTimeSolution, bc,elements,timeStep)

tol = 1e-2;
iter = 0;
maxIter = 50;
error = 100;
errorT = 100;
maxLineSearchIter = 5;
while error > tol && iter < maxIter && errorT > 0.3

    [J, R] = residualCalculation(prevSolution, prevTimeSolution, bc, elements, timeStep);
    dX = J\R;
    error = norm(R,2);
    
    iterLineSearch = 0;
    alpha = 1.0;
    while iterLineSearch < maxLineSearchIter 
        solution = prevSolution + alpha*dX;
        % function  to evaluate residual with new solution
        [~, R] = residualCalculation(solution, prevTimeSolution, bc, elements, timeStep);
        errorNew = norm(R,2);

        if errorNew > error
            alpha = alpha * 0.5;
            iterLineSearch = iterLineSearch + 1;
        else
            break;
        end
    
    end
    
    errorT = norm(solution - prevSolution,2);
    
    prevSolution = solution;
    iter = iter + 1;
    
end