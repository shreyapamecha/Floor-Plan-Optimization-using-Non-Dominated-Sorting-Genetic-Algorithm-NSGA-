# Floor-Plan-Optimization-using-Non-Dominated-Sorting-Genetic-Algorithm-NSGA-

The floor Plan Optimization problem has two conflicting objectives: Area of the Laboratory (Maximizing
function) and Cost of designing the workspace (Minimizing function). To solve this, we need to use a
multi-objective decision-making method: Non-Dominated Sorting Genetic Algorithm.

Number of particles chosen = 10
Maximum Iterations = 30
The new Chemistry Workspace has 3 rooms - laboratory, equipment room, and sitting/meeting room.

According to the constraints given in the question,
● Laboratory: 10 <= x <= 15, 10 <= 25-y <= 18
● Sitting/Meeting room: length=25, 10 <= x <= 15
● Equipment Room: 10 <= 25-x <=15, 7 <= y <=15

Particles Size (10 x 8) - where 10 is the number of particles and 8 is the number of columns containing
binary digits to represent 2 integer numbers from 0 to 15 for x & y. (Could have increased the number of
binary digits to represent decimal numbers for x & y).

Methodology:
1. Initialize all the 10 particles (x,y) in the ‘pop’ matrix.
2. Decipher the integer values for x and y, then the functional values (area, cost), and put them in
the matrix in columns 9, 10, 11, and 12 respectively.
3. Determine their ranks using Non-Dominated Sorting and put it as column ‘13’.
4. Sort the ‘pop’ matrix with respect to the 13th column containing the ranks.
5. Now, we need to select 8 parents (Crossover Probability = 0.8) for the Crossover operation.
Before that, we need to find whether we require to sort the particles within a rank.
6. If yes, then, use the crowding distance method to sort the particles within a rank and put it in the
main matrix.
7. Select the top 8 particles randomly for crossover operation and generate 8 children.
8. Determine the functional values for all these 8 children, assign all the 18 particles proper ranks,
and sort the matrix with respect to these ranks.
9. Now, select the top 10 particles after employing crowding distance to sort the particles within a
rank.
10. Repeat this process.

Video Link: https://youtu.be/VeRFErNBmdE
The video shared shows different Pareto Fronts or set of optimal solutions with increasing iterations in
multi-objective optimization problems stands for a set of solutions that are non-dominated to each other
but are superior to the rest of the solutions in the search space.
In the end, we obtain a straight line that conveys the information of the sheer trade-off between the two
objective functions.
