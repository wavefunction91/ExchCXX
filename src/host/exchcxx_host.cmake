file(GLOB EXCHCXX_HOST_SOURCES "cpu/*.cpp")
target_sources( exchcxx PRIVATE ${EXCHCXX_HOST_SOURCES} )
