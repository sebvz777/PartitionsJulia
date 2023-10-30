@test_set "PartitionProperties test" begin
    @test_set "pair" begin
        @test is_pair(Partition([1, 2, 2, 1], [3, 3]))
        @test !is_pair(Partition([1, 2, 2, 1], [4, 3]))
        @test !is_pair(Partition([1, 1, 1], [3, 3, 3, 3]))
        @test is_pair(Partition([2, 1], [1, 2]))
    end
    @test_set "balanced" begin
        @test is_balanced(Partition([1, 2, 2, 1], [3, 3]))
        @test is_balanced(Partition([1, 2, 3], [3, 2, 1]))
        @test !is_balanced(Partition([1, 1, 3], [3, 3]))
    end
    @test_set "non-crossing" begin
        @test is_noncrossing(Partition([1, 2, 3, 4, 5, 5, 4, 3, 2, 1], [1, 6, 6, 1]))
        @test !is_noncrossing(Partition([1, 2, 3, 4, 5, 5, 4, 3, 2, 1], [1, 6, 6, 2]))
    end
end