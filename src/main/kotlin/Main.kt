@file:JvmName("Main")

import org.apache.commons.math3.special.Erf
import kotlin.math.log10
import kotlin.math.log2
import kotlin.math.sqrt
import kotlin.random.Random

data class Matrix(val source: List<List<Int>>) {
    fun height() = source.size
    fun width() = source.first().size

    constructor(vararg rows: List<Int>) : this(rows.toList())

    override fun toString(): String =
        source.joinToString(prefix = "[", postfix = "]", separator = "\n ") { Vector(it).toString() }
}

fun Matrix.transpose(): Matrix {
    return Matrix(
        (0 until width()).map { j ->
            (0 until height()).map { i ->
                source[i][j]
            }
        }
    )
}

data class Vector(val source: List<Int>) {
    fun length() = source.size

    constructor(vararg args: Int) : this(args.toList())

    override fun toString(): String = source.joinToString(prefix = "[", postfix = "]") { it.toString() }
}

operator fun Vector.times(matrix: Matrix): Vector {
    if (matrix.height() != length()) error("Invalid operation vector${length()} * matrix(${matrix.height()} x ${matrix.width()})")
    return Vector(
        (0 until matrix.width()).map { i ->
            (0 until length()).sumBy { j ->
                source[j] * matrix.source[j][i]
            } % 2
        }
    )
}

fun allVectorsSize(size: Int): List<Vector> {
    if (size == 0) return emptyList()
    if (size == 1) return listOf(Vector(listOf(0)), Vector(listOf(1)))
    return allVectorsSize(size - 1).flatMap {
        listOf(Vector(it.source + 0), Vector(it.source + 1))
    }
}

val G = Matrix(
    listOf(1, 0, 0, 0, 0, 0, 0, 0, 1, 1),
    listOf(0, 1, 0, 0, 0, 0, 0, 1, 1, 1),
    listOf(0, 0, 1, 0, 0, 0, 1, 1, 0, 1),
    listOf(0, 0, 0, 1, 0, 0, 1, 1, 1, 0),
    listOf(0, 0, 0, 0, 1, 0, 1, 1, 1, 1),
    listOf(0, 0, 0, 0, 0, 1, 1, 1, 0, 0)
)

val HT = Matrix(
    listOf(0, 1, 1, 1, 0, 0, 1, 1, 1, 1),
    listOf(0, 1, 0, 0, 1, 1, 0, 1, 1, 1),
    listOf(1, 0, 1, 0, 0, 1, 0, 1, 1, 0),
    listOf(0, 0, 0, 0, 1, 1, 1, 0, 1, 1)
).transpose()

const val SPEED = 0.6
const val DEFAULT_COUNT = 1000

val CODES = allVectorsSize(6).associateBy { it * G }

infix fun Vector.hammingTo(other: Vector): Int {
    if (length() != other.length()) error("Error 15") // lol
    return source.zip(other.source).sumBy { (a, b) -> 0.takeIf { a == b } ?: 1 }
}

infix fun Vector.closestFrom(vectors: Collection<Vector>): Vector {
    return vectors.minBy { it hammingTo this }!!
}

operator fun Matrix.times(other: Matrix): Matrix {
    if (width() != other.height()) error("Error 16")
    return Matrix(
        (0 until height()).map { i ->
            (0 until other.width()).map { j ->
                (0 until width()).sumBy { s ->
                    source[i][s] * other.source[s][j]
                } % 2
            }
        }
    )
}

infix fun Vector.withError(probability: Double): Vector {
    return Vector(source.map { it.takeIf { Random.nextDouble() >= probability } ?: (it + 1) % 2 })
}

fun errorProbability(modulationLevels: Int, noise: Double, speed: Double): Double {
    val error = bitErrorRate(noise, modulationLevels, speed)
    return allVectorsSize(6).map { source ->
        (0 until DEFAULT_COUNT).map {
            val received = (source * G) withError error
            val decoded = CODES[(received closestFrom CODES.keys)]
            if (decoded != source) {
                1
            } else {
                0
            }
        }.average()
    }.average()
}

fun rawErrorProbability(modulationLevels: Int, noise: Double): Double {
    val error = bitErrorRate(noise, modulationLevels, 1.0)
    return allVectorsSize(6).map { source ->
        (0 until DEFAULT_COUNT).map {
            if ((source withError error) != source) {
                1
            } else {
                0
            }
        }.average()
    }.average()
}

// https://old.telesputnik.ru/archive/pdf/181/70.pdf
// approximation, error ~ 0.2db
fun bitErrorRate(noise: Double, modulationLevels: Int, speed: Double): Double {
    val energy = noise - 10 * log10(log2(modulationLevels.toDouble() / speed))
    return 2 * (1 - 1 / sqrt(modulationLevels.toDouble())) *
            Erf.erfc(sqrt(3 * log2(modulationLevels.toDouble()) * energy / (2 * (modulationLevels - 1)))) / energy
}

fun binarySearch(left: Double, right: Double, choose: (Double) -> Boolean): Double {
    var l = left
    var r = right
    for (i in 1..32) {
        val middle = ((l + r) / 2).coerceIn(l, r)
        if (choose(middle)) {
            r = middle
        } else {
            l = middle
        }
    }
    return l
}

fun main() {
    listOf(2, 4, 8, 16, 64).forEach br@{ levels ->
        println()
        (1..100).step(1).forEach { noise ->
            if (!bitErrorRate(noise.toDouble(), levels, SPEED).isNaN()) {
                val error = errorProbability(levels, noise.toDouble(), SPEED)
                println("[${levels}QAM] signal/noise=${noise}dB, \terror=$error")
                if (error < 1e-10) return@br
            }
        }
    }

    listOf(2, 4, 8, 16, 64).forEach { levels ->
        val raw = binarySearch(0.0, 100.0) {
            if (it.isNaN()) return@binarySearch false
            rawErrorProbability(levels, it) < 1e-5
        }
        val cool = binarySearch(0.0, 100.0) {
            if (it.isNaN()) return@binarySearch false
            errorProbability(levels, it, SPEED) < 1e-5
        }
        println("[${levels}QAM] p=${1e-5}, efficiency(S'/N' - S/N)=${raw - cool}")
    }

}