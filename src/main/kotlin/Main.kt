@file:JvmName("Main")

import org.apache.commons.math3.distribution.NormalDistribution
import kotlin.math.pow
import kotlin.math.sqrt

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

data class VectorDouble(val source: List<Double>)

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
val ALL_VECTORS = allVectorsSize(6)

infix fun Vector.closestFromByHamming(vectors: Collection<Vector>): Vector {
    infix fun Vector.hammingTo(other: Vector): Int {
        if (length() != other.length()) error("Error 15") // lol
        return source.zip(other.source).sumBy { (a, b) -> 0.takeIf { a == b } ?: 1 }
    }
    return vectors.minBy { it hammingTo this }!!
}

infix fun VectorDouble.closestFrom(vectors: Collection<Vector>): Vector {
    return vectors.maxBy { this * it }!!
}

operator fun VectorDouble.times(vector: Vector): Double {
    return source.zip(vector.source).sumByDouble { (a, b) -> a * b }
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

infix fun Vector.withNoise(distribution: NormalDistribution): VectorDouble {
    return VectorDouble(source.map { 2 * it - 1 + distribution.sample() })
}

infix fun Vector.withHardNoise(distribution: NormalDistribution): Vector {
    return Vector(source.map { if (2 * it - 1 + distribution.sample() < 0) 0 else 1 })
}

fun hardErrorProbability(noise: NormalDistribution): Double {
    return ALL_VECTORS.map { source ->
        (0 until DEFAULT_COUNT).map {
            val received = (source * G) withHardNoise noise
            val decoded = CODES[received closestFromByHamming CODES.keys]
            if (decoded != source) {
                1
            } else {
                0
            }
        }.average()
    }.average()
}

fun softErrorProbability(noise: NormalDistribution): Double {
    return ALL_VECTORS.map { source ->
        (0 until DEFAULT_COUNT).map {
            val received = (source * G) withNoise noise
            val decoded = CODES[(received closestFrom CODES.keys)]
            if (decoded != source) {
                1
            } else {
                0
            }
        }.average()
    }.average()
}

fun rawErrorProbability(noise: NormalDistribution): Double {
    return ALL_VECTORS.map { source ->
        (0 until DEFAULT_COUNT).map {
            if (((source withNoise noise) closestFrom ALL_VECTORS) != source) {
                1
            } else {
                0
            }
        }.average()
    }.average()
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

fun sigma(noise: Double): Double = sqrt(0.5 * 10.0.pow(-noise / 10.0) / SPEED)

fun main() {
    (1..30).step(1).forEach { noise ->
        val error = softErrorProbability(NormalDistribution(0.0, sigma(noise.toDouble())))
        println("[SOFT] Eb/N0=${noise}dB, \terror=$error")
    }

    (1..30).step(1).forEach { noise ->
        val error = hardErrorProbability(NormalDistribution(0.0, sigma(noise.toDouble())))
        println("[HARD] Eb/N0=${noise}dB, \terror=$error")
    }

    println("P=${1e-5}")

    val raw = binarySearch(0.0, 100.0) {
        rawErrorProbability(NormalDistribution(0.0, sigma(it))) < 1e-5
    }
    println("[NO DECODING]: Eb/N0 = $raw")
    val hard = binarySearch(0.0, 100.0) {
        hardErrorProbability(NormalDistribution(0.0, sigma(it))) < 1e-5
    }
    println("[HARD]: Eb/N0 = $hard, efficiency = ${raw - hard}")
    val cool = binarySearch(0.0, 100.0) {
        softErrorProbability(NormalDistribution(0.0, sigma(it))) < 1e-5
    }
    println("[SOFT]: Eb/N0 = $raw, efficiency = ${raw - cool}")

}