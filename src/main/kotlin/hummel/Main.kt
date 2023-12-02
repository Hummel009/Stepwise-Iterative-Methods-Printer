package hummel

import java.text.DecimalFormat
import java.util.*
import kotlin.math.abs

fun main() {
	var method: IterationsMethod
	val m = Array(10) { DoubleArray(11) }
	for (i in 0..9) {
		for (j in 0..9) {
			if (i == j) {
				m[i][j] = 20.0
			} else {
				m[i][j] = 1.0
			}
		}
		m[i][10] = (19 * (i + 1) + 171).toDouble()
	}
	method = IterationsMethod(m, true)
	println("Comparison of Jacobi and Gauss-Seidel methods for System of linear equations")
	println()
	println("Source matrix:")
	println()
	method.print()
	println()
	println("Jacobi Iterations:")
	println()
	method.solve()
	println()
	println("Seidel-Gauss Iterations:")
	println()
	method = IterationsMethod(m, false)
	method.solve()
}

class IterationsMethod(private var m: Array<DoubleArray>, private var isJacobi: Boolean) {
	private val maxIterations = 10000

	@Suppress("LoopToCallChain")
	fun print() {
		val n = m.size
		for (doubles in m) {
			for (j in 0 until n + 1) {
				val result = DecimalFormat("#.#####").format(doubles[j])
				print("$result ")
			}
			println()
		}
	}

	@Suppress("LoopToCallChain")
	fun solve() {
		var iterations: Int = if (isJacobi) 1 else 0

		val n = m.size
		val epsilon = 1e-4
		val x = DoubleArray(n)
		var p = DoubleArray(n)
		Arrays.fill(x, 0.0)
		if (isJacobi) {
			Arrays.fill(p, 0.0)
		}
		while (true) {
			for (i in 0 until n) {
				var sum = m[i][n]
				for (j in 0 until n) {
					if (j != i) {
						sum -= if (isJacobi) {
							m[i][j] * p[j]
						} else {
							m[i][j] * x[j]
						}
					}
				}
				x[i] = 1 / m[i][i] * sum
			}
			print("K = $iterations; X: ")
			for (i in 0 until n) {
				val result = DecimalFormat("#.#####").format(x[i])
				print("$result; ")
			}
			println()
			if (iterations != 0) {
				val result = DecimalFormat("#.#####").format(abs(x[0] - p[0]))
				println("||X(" + iterations + ")-X(" + (iterations - 1) + ")|| = " + result)
			}
			println()
			iterations++
			if (iterations == 1) {
				continue
			}
			var stop = true
			for (i in 0 until n) {
				if (abs(x[i] - p[i]) > epsilon) {
					stop = false
					break
				}
			}
			if (stop || iterations == maxIterations) {
				break
			}
			p = x.clone()
		}
	}
}