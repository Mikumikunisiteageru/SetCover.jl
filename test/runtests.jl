using SetCover
using Test
occurrences = [1 2; 1 3; 1 4; 1 5; 2 1; 2 2; 2 3; 3 4; 3 5; 3 6]
areas = [1 1.3; 2 1.0; 3 1.0]
scp = SCP(occurrences, areas)
@test complementary(scp).solution == [1, 1, 1]
@test        greedy(scp).solution == [1, 1, 1]
@test    lagrangian(scp).solution == [0, 1, 1]
@test    ilp_solver(scp).solution == [0, 1, 1]
