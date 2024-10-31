struct CheybyshevT{N} end
(ch::CheybyshevT{0})(x::T) where T = 1.0
(ch::CheybyshevT{1})(x::T) where T = x
(ch::CheybyshevT{2})(x::T) where T = 2x^2 - 1
(ch::CheybyshevT{3})(x::T) where T = 4x^3 - 3x
(ch::CheybyshevT{4})(x::T) where T = 8x^4 - 8x^2 + 1
(ch::CheybyshevT{5})(x::T) where T = 16x^5 - 20x^3 + 5x
