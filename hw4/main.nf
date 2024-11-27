
// Задаём параметры через params
params.outdir = "./test_output"  // Путь для вывода по умолчанию
params.threads = 8               // Количество потоков по умолчанию

params.spades = [
    threads: 4,                  // Параметры для SPADES
    memory: "8G"
]

params.quast = [
    reference: ""  // Путь к референсному файлу
]

params.prokka = [
    genus: "Escherichia"          // Род по умолчанию
]

params.abricate = [
    database: "ncbi"              // База данных по умолчанию для Abricate
]

import groovy.json.JsonSlurper

// Функция для загрузки JSON и обновления параметров
def loadParamsFromJson(jsonPath) {
    def jsonFile = new File(jsonPath)
    def json = new JsonSlurper().parseText(jsonFile.text)
    
    // Обновляем параметры через params
    params.outdir = json.global_params.outdir ?: params.outdir
    params.threads = json.global_params.threads ?: params.threads
    
    params.spades.threads = json.spades?.threads ?: params.spades.threads
    params.spades.memory = json.spades?.memory ?: params.spades.memory
    
    params.quast.reference = json.quast?.reference ?: params.quast.reference
    
    params.prokka.genus = json.prokka?.genus ?: params.prokka.genus
    
    params.abricate.database = json.abricate?.database ?: params.abricate.database
}

// Загружаем параметры из JSON файла
loadParamsFromJson('./params.json')

// Процесс для FastQC
process fastqc_n {
    tag { sample_name }
    publishDir("${params.outdir}/fastqc", mode: "copy")
    input:
    tuple val(sample_name), path(read1), path(read2), val(assembly)

    output:
    path "${sample_name}_fastqc/", emit: fastq_n_o

    script:
    """
    mkdir -p ${sample_name}_fastqc
    fastqc -o ${sample_name}_fastqc ${read1} ${read2}
    """
}

// Процесс для SPADES
process spades_n {
    tag {sample_name}

    publishDir("${params.outdir}/spades", mode: "copy")
    input:
    tuple val(sample_name), path(read1), path(read2), val(assembly)
    
    output:
    path "${sample_name}_spades/", emit: spades_n_o

    script:
    """
    mkdir -p ${sample_name}_spades
    spades.py -1 ${read1} -2 ${read2} -o ${sample_name}_spades --threads ${params.spades.threads} --memory ${params.spades.memory.replaceAll("[^\\d]", "")}

    """
}
process quast_n{
    tag {sample_name}

    publishDir ("./test_output/quast", mode: "copy")
    input:
    path assembly_path
    tuple val(sample_name), path(read1), path(read2), val(assembly)
    
    output:
    path "${sample_name}_quast/", emit: quast_n_o

    script:
    """
    mkdir -p ${sample_name}_quast
    quast.py ${assembly_path}/scaffolds.fasta -o ${sample_name}_quast
    """
}

process prokka_n {
    tag { sample_name }

    publishDir("./test_output/prokka", mode: "copy")
    
    input:
    path assembly_path // путь к сборке генома
    tuple val(sample_name), path(read1), path(read2), val(assembly)

    output:
    path "${sample_name}_prokka/", emit: prokka_n_o

    script:
    """
    mkdir -p ${sample_name}_prokka
    prokka --outdir ${sample_name}_prokka --prefix ${sample_name} --force ${assembly_path}/scaffolds.fasta 
    """
}
process abricate_n {
    tag { sample_name }

    publishDir("./test_output/abricate", mode: "copy")
    
    input:
    path assembly_path // путь к сборке генома
    tuple val(sample_name), path(read1), path(read2), val(assembly)

    output:
    path "${sample_name}_abricate/", emit: abricate_n_o

    script:
    """
    mkdir -p ${sample_name}_abricate
    abricate --db ncbi ${assembly_path}/scaffolds.fasta > ${sample_name}_abricate/abricate_results.txt
    """
}



process quast{
    tag {sample_name}

    publishDir ("./test_output/quast", mode: "copy")
    input:
    tuple val(sample_name), path(assembly)
    
    output:
    path "${sample_name}_quast/", emit: quast_n_o

    script:
    """
    mkdir -p ${sample_name}_quast
    quast.py ${assembly} -o ${sample_name}_quast
    """
}

process prokka {
    tag { sample_name }

    publishDir("./test_output/prokka", mode: "copy")
    
    input:
    tuple val(sample_name), path(assembly)

    output:
    path "${sample_name}_prokka/", emit: prokka_n_o

    script:
    """
    mkdir -p ${sample_name}_prokka
    prokka --outdir ${sample_name}_prokka --prefix ${sample_name} --force ${assembly}
    """
}
process abricate {
    tag { sample_name }

    publishDir("./test_output/abricate", mode: "copy")
    
    input:
    tuple val(sample_name), path(assembly)

    output:
    path "${sample_name}_abricate/", emit: abricate_n_o

    script:
    """
    mkdir -p ${sample_name}_abricate
    abricate --db ncbi ${assembly} > ${sample_name}_abricate/abricate_results.txt
    """
}


// Чтение CSV и выполнение FastQC
workflow {

    def samples = Channel.fromPath('./samples.csv')  // Путь к вашему CSV
        .splitCsv(header: true)  // Разбиваем на строки с учетом заголовка
        .map { row -> 
            def sample_id = row.sample_id
            def read1 = file(row.read_1)
            def read2 = file(row.read_2)
            def assembly = row.assembly // Считываем поле assembly, если оно есть
            return [sample_id, [read1, read2, assembly]]  // Возвращаем пару ключ-значение
        }

  
    sample_list = samples.map{ sample_name, reads -> 
    def combined = [sample_name] + reads
    return combined
    }


    sample_list.branch { item ->
    with_scaffolds: item[-1] != ""
    without_scaffolds: item[-1] == ""
}.set { result }


// result.with_scaffolds.view { println "With scaffolds: $it" }
def assembly_channel = result.with_scaffolds.map { item -> 
    return [item[0], file(item[3])]  // item[0] — sample_name, item[3] — сборка
}
// assembly_channel.view()
// result.without_scaffolds.view { println "Without scaffolds: $it" }


fastqc_n(result.without_scaffolds)
spades_n(result.without_scaffolds)
spades_n.out.spades_n_o.view()
quast_n(spades_n.out.spades_n_o, result.without_scaffolds )
prokka_n(spades_n.out.spades_n_o, result.without_scaffolds)
abricate_n(spades_n.out.spades_n_o, result.without_scaffolds)

quast(assembly_channel)
prokka(assembly_channel)
abricate(assembly_channel)

// result.without_scaffolds.view { items ->
//     if (items.size() > 0) {
//         println "Без скаффолда есть"
//     }
// }
// result.with_scaffolds.view { items ->
//     if (items.size() > 0) {
//         println "С скаффолдом есть"
//     }
// }

}

