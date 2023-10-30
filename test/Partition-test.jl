@testset "Partition.jl test" begin
    @testset "Partition constructor" begin
        @test Partition([2, 3], [8, 9]).upper_points == [1, 2] && Partition([2, 3], [8, 9]).lower_points == [3, 4]
    end
    @testset "Operations" begin
        @test involution(Partition([2, 3], [4, 4])) == Partition([1, 1], [2, 3])
        @test tensor_product(Partition([1, 2, 3], [2, 1, 2]), Partition([3, 3], [])) == Partition([1, 2, 3, 4, 4], [2, 1, 2])
        @test_throws composition(Partition([1, 2, 3], [3]), Partition([1, 2], [1, 2]))
        @test composition(Partition([1, 2, 2, 2], [2, 3, 4]), Partition([1, 2], [3, 3, 1, 2])) == Partition([1, 1], [1, 2, 3])
        @test_throws rotation([], [1], true, true)
        @test rotation([1, 2], [2, 1], false, false) == Partition([1, 2, 1], [2])
        @test vertical_reflection(Partition([2, 3, 2, 2], [2, 3])) == Partition([1, 1, 2, 1], [2, 1])
    end
end