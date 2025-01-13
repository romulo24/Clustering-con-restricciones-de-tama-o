# Librerías necesarias
library(cluster)
library(proxy)
library(mlr3oml)
library(mlr3)
library(pryr)
library(dplyr) 
library(aricode)
library(ggplot2)
library(corrplot)
library(clValid)
library(RColorBrewer)


# Función para preparar datos con validación de variable dependiente
prepare_data <- function(dataset) {
  # Asegurarse de que el dataset es un data.frame
  dataset <- as.data.frame(dataset)
  
  # Verificar si existe la columna "class"
  if ("class" %in% colnames(dataset)) {
    # Separar la variable "class" y el resto de las columnas
    y <- dataset$class
    X <- dataset[, setdiff(colnames(dataset), "class")]
  } 
  # Validar si la última columna es categórica o entera
  else if (is.factor(dataset[, ncol(dataset)]) || is.integer(dataset[, ncol(dataset)])) {
    X <- dataset[, -ncol(dataset)]
    y <- dataset[, ncol(dataset)]
  } 
  # Si ninguna condición se cumple, devolver NULL para ignorar este dataset
  else {
    return(NULL)
  }
  
  # Devolver una lista con las variables X e y
  list(X = X, y = y)
}

# Función para calcular la matriz de distancia
calculate_distance_matrix <- function(X) {
  D <- proxy::dist(as.matrix(X), method = "euclidean")
  as.matrix(D)
}

# Función para ajustar cardinalidad
adjust_cardinality <- function(cluster_assignment, X, centroids, target_cardinality) {
  max_iterations <- 1000 # Evitar bucles infinitos
  iteration <- 0
  
  for (j in 1:length(target_cardinality)) {
    cluster_size <- sum(cluster_assignment == j, na.rm = TRUE)
    
    while (cluster_size > target_cardinality[j] && iteration < max_iterations) {
      # Incrementar el contador de iteraciones
      iteration <- iteration + 1
      
      # Obtener índices de elementos en el cluster excedente
      idx <- which(cluster_assignment == j)
      element <- idx[1]
      
      # Calcular distancias del elemento a todos los centroides
      distances <- apply(centroids, 1, function(centroid) sum((X[element, ] - centroid)^2))
      
      # Identificar clusters con espacio disponible
      available_clusters <- which(sapply(1:length(target_cardinality), function(c) {
        sum(cluster_assignment == c, na.rm = TRUE) < target_cardinality[c]
      }))
      
      # Filtrar los clusters alternativos válidos
      valid_clusters <- available_clusters[available_clusters != j]
      
      # Si no hay clusters válidos, marcar como "sin asignar" (valor especial)
      if (length(valid_clusters) == 0) {
        cluster_assignment[element] <- -1
        cluster_size <- sum(cluster_assignment == j, na.rm = TRUE)
        next
      }
      
      # Ordenar los clusters válidos por cercanía
      sorted_clusters <- valid_clusters[order(distances[valid_clusters])]
      
      # Reasignar el elemento al cluster más cercano con espacio
      cluster_assignment[element] <- sorted_clusters[1]
      
      # Actualizar el tamaño del cluster actual
      cluster_size <- sum(cluster_assignment == j, na.rm = TRUE)
    }
  }
  
  # Si hay elementos sin asignar, reubicarlos al cluster con menor tamaño
  if (any(cluster_assignment == -1)) {
    unassigned_idx <- which(cluster_assignment == -1)
    for (element in unassigned_idx) {
      # Encontrar el cluster con espacio disponible y menor tamaño
      smallest_valid_cluster <- which.min(sapply(1:length(target_cardinality), function(c) {
        if (sum(cluster_assignment == c, na.rm = TRUE) < target_cardinality[c]) {
          sum(cluster_assignment == c, na.rm = TRUE)
        } else {
          Inf
        }
      }))
      
      # Asignar el elemento al cluster con menor tamaño
      cluster_assignment[element] <- smallest_valid_cluster
    }
  }
  
  # Mensaje de advertencia si se alcanzó el límite de iteraciones
  if (iteration >= max_iterations) {
    warning("El proceso de ajuste alcanzó el máximo de iteraciones.")
  }
  
  return(cluster_assignment)
}

# Función para generar solución inicial
generate_initial_solution <- function(X, target_cardinality, seed = 45) {
  set.seed(seed)
  km <- kmeans(X, centers = length(target_cardinality))
  cluster_assignment <- km$cluster
  centroids <- calculate_centroids(X, cluster_assignment, length(target_cardinality))
  adjust_cardinality(cluster_assignment, X, centroids, target_cardinality)
}

# Función para calcular centroides
calculate_centroids <- function(X, cluster_assignment, k) {
  centroids <- matrix(0, nrow = k, ncol = ncol(X))
  for (j in 1:k) {
    centroids[j, ] <- colMeans(X[cluster_assignment == j, , drop = FALSE])
  }
  return(centroids)
}

# Función para evaluar una solución
evaluate_solution <- function(cluster_assignment, X, target_cardinality, penalty_weight = 10) {
  ss <- silhouette(cluster_assignment, dist(X))
  penalty <- 0
  for (j in 1:length(target_cardinality)) {
    cardinality_diff <- abs(sum(cluster_assignment == j) - target_cardinality[j])
    penalty <- penalty + cardinality_diff * penalty_weight
  }
  return(mean(ss[, "sil_width"]) - penalty)
}

# Algoritmo de murciélagos
run_bat_algorithm <- function(X, y, target_cardinality, n_bats = 30, max_iterations = 20,
                              f_min = 0, f_max = 2, loudness = 0.5, pulse_rate = 0.5,
                              alpha = 0.9, gamma = 0.9) {
  set.seed(1521) #set.seed(1807)
  seeds <- sample(1:10000, n_bats, replace = FALSE)
  bats <- lapply(1:n_bats, function(i) {
    list(
      position = generate_initial_solution(X, target_cardinality, seed = seeds[i]),
      velocity = rep(0, nrow(X)),
      frequency = runif(1, f_min, f_max),
      loudness = loudness,
      pulse_rate = pulse_rate,
      seed = seeds[i]
    )
  })
  best_solution <- bats[[1]]$position
  best_score <- evaluate_solution(best_solution, X, target_cardinality)
  best_seed <- bats[[1]]$seed
  for (iteration in 1:max_iterations) {
    for (i in 1:n_bats) {
      bat <- bats[[i]]
      bat$frequency <- runif(1, f_min, f_max)
      bat$velocity <- bat$velocity + (bat$position - best_solution) * bat$frequency
      new_position <- round(bat$position + bat$velocity)
      new_position <- pmin(pmax(new_position, 1), length(target_cardinality))
      if (runif(1) > bat$pulse_rate) {
        new_position <- sample(1:length(target_cardinality), nrow(X), replace = TRUE)
      }
      for (j in 1:length(target_cardinality)) {
        while (sum(new_position == j) > target_cardinality[j]) {
          idx <- which(new_position == j)
          new_position[sample(idx, 1)] <- sample(1:length(target_cardinality), 1)
        }
      }
      new_score <- evaluate_solution(new_position, X, target_cardinality)
      if (new_score > best_score && runif(1) < bat$loudness) {
        bats[[i]]$position <- new_position
        bats[[i]]$loudness <- alpha * bat$loudness
        bats[[i]]$pulse_rate <- pulse_rate * (1 - exp(-gamma * iteration))
        if (new_score > best_score) {
          best_solution <- new_position
          best_score <- new_score
          best_seed <- bats[[i]]$seed
        }
      }
    }
  }
  list(best_solution = best_solution, best_score = best_score, best_seed = best_seed, seeds = seeds)
}


print_results <- function(results, y, X, D, target_cardinality, dataset_name) {
  best_solution <- results$best_solution
  best_score <- results$best_score
  best_seed <- results$best_seed
  seeds <- results$seeds
  num_instancias <- nrow(X)  
  num_variables <- ncol(X) + 1
  
  # Calcular ARI, AMI y NMI
  ARI_value <- ARI(y, best_solution)
  AMI_value <- AMI(y, best_solution)
  NMI_value <- NMI(y, best_solution)
  
  # Calcular el coeficiente de silueta para la solución final
  silhouette_values <- silhouette(x = best_solution, dist = as.dist(D))
  mean_silhouette <- mean(silhouette_values[, "sil_width"])
  
  # Convertir los valores de silueta en un data frame
  silhouette_df <- as.data.frame(silhouette_values[, c("cluster", "sil_width")])
  
  # Guardar como archivo CSV
  write.csv(silhouette_df, "silhouette_results.csv", row.names = FALSE)
  
  cat("Los resultados del coeficiente de silueta se han guardado en 'silhouette_results.csv'.")
  
  #Clacular el indice de dunn
  dunn_index <- dunn(distance = as.dist(D), clusters = best_solution)

  # Contar el número de clústeres
  num_clusters <- length(unique(best_solution))
  class_dist = as.integer(table(best_solution))
  # Almacenar los resultados en el data frame global
  global_results <<- rbind(global_results, data.frame(
    name = dataset_name,
    Best_Seed = best_seed,
    ARI = ARI_value,
    AMI = AMI_value,
    NMI = NMI_value,
    Mean_Silhouette = mean_silhouette,
    Dunn_index = dunn_index,
    Clusters = num_clusters,
    number_features = num_variables,
    number_instances= num_instancias,
    cardinality_BAT = I(list(class_dist)),
    cardinality_REAL = I(list(target_cardinality))
    
  ))
  
  # Crear una paleta de colores pastel
  colores_pastel <- brewer.pal(n = max(best_solution), name = "Pastel1")
  
  # Graficar los valores de silueta
  plot(silhouette_values, col = colores_pastel, border = NA, cex.names = 0.7)
  abline(v = mean_silhouette, lty = 2, col = "red")
  text(mean_silhouette, max(silhouette_values[, "sil_width"]),
       labels = paste("Media =", round(mean_silhouette, 3)),
       pos = 2, col = "red", cex = 0.8)
  
  # Mostrar información en consola
  cat("Promedio del coeficiente de silueta:", mean_silhouette, "\n")
  cat("Indice de Dunn:", dunn_index, "\n")
  cat("Semillas usadas para cada murciélago:\n")
  print(seeds)
  cat("\nLa semilla del murciélago con la mejor solución fue:", best_seed, "\n")
  
  cat("\nAsignación de clusters óptima (Algoritmo de Murciélagos):\n")
  print(table(best_solution))
  cat("Número de clusters:", num_clusters, "\n")
  
  cat("\nAdjusted Rand Index (ARI):", ARI_value, "\n")
  cat("Adjusted Mutual Information (AMI):", AMI_value, "\n")
  cat("Normalized Mutual Information (NMI):", NMI_value, "\n")
}


# Función principal para correr todo
run_clustering <- function(dataset, target_cardinality, dataset_name) {
  data <- prepare_data(dataset)
  X <- data$X
  y <- data$y
  D <- calculate_distance_matrix(X)
  results <- run_bat_algorithm(X, y, target_cardinality)
  print_results(results, y, X, D, target_cardinality, dataset_name)
  
}



# Función para filtrar datasets únicos
filtrar_datasets_unicos <- function(dataset) {
  dataset %>%
    distinct(NumberOfClasses, NumberOfInstances, MajorityClassSize, MinorityClassSize, .keep_all = TRUE) %>%
    mutate(name_lower = tolower(name)) %>%
    distinct(name_lower, .keep_all = TRUE) %>%
    select(-name_lower) %>%
    filter(!is.na(NumberOfSymbolicFeatures), NumberOfSymbolicFeatures <= 1) %>%
    filter(!is.na(NumberOfMissingValues), NumberOfMissingValues == 0)
}

# Función para cargar datasets con manejo de errores
cargar_dataset <- function(id) {
  tryCatch({
    odata <- OMLData$new(id = id)
    dataset <- odata$data
    
    if (any(duplicated(names(dataset)))) {
      warning(paste("Dataset ID", id, ": Columnas duplicadas encontradas. Renombrando columnas..."))
      names(dataset) <- make.unique(names(dataset))
    }
    
    message(paste("Dataset ID", id, "descargado con éxito."))
    return(dataset)
  }, error = function(e) {
    warning(paste("Error al cargar el dataset ID", id, ":", e$message))
    return(NULL)
  })
}

# Función para obtener distribuciones de clases y validar cantidad de clases
get_class_distributions <- function(id, expected_classes) {
  dataset <- cargar_dataset(id)
  
  # Si el dataset no pudo descargarse, devolver valores nulos
  if (is.null(dataset)) {
    return(list(text = "Error", vector = NA, data = NULL))
  }
  
  # Determinar la variable dependiente (última columna o "class")
  if ("class" %in% colnames(dataset)) {
    target <- dataset$class
  } else {
    target <- dataset[[ncol(dataset)]]  # Última columna
  }
  
  # Obtener la distribución de las clases
  class_dist <- table(target)
  
  # Validar cantidad de clases
  num_classes_found <- length(class_dist)
  if (num_classes_found != expected_classes) {
    warning(paste("Dataset ID", id, "tiene una discrepancia en el número de clases:",
                  "esperado =", expected_classes, ", encontrado =", num_classes_found))
    return(list(text = "Error", vector = NA, data = NULL))
  }
  
  # Formato de texto
  class_dist_text <- paste(names(class_dist), as.integer(class_dist), sep = ":", collapse = "; ")
  
  # Formato de vector numérico
  class_dist_vector <- as.integer(class_dist)
  
  # Devolver lista con la distribución de clases y el dataset
  return(list(text = class_dist_text, vector = class_dist_vector, data = dataset))
}

# Obtener lista de datasets
odatasets <- list_oml_data(
  number_features = c(2, 15),
  number_instances = c(100, 5000),
  number_classes = c(2, 15)
)

# Filtrar datasets únicos
odatasets_unique <- filtrar_datasets_unicos(odatasets)

# Medir el tiempo de descarga
start_time <- Sys.time()

# Recorrer cada data_id y validar contra la cantidad esperada de clases
class_distributions <- mapply(
  FUN = get_class_distributions,
  id = odatasets_unique$data_id,
  expected_classes = odatasets_unique$NumberOfClasses,
  SIMPLIFY = FALSE
)

end_time <- Sys.time()

execution_time <- end_time - start_time
print(paste("Tiempo total de ejecución:", execution_time))

# Filtrar datasets válidos
valid_datasets <- Filter(function(x) !is.null(x$data), class_distributions)

# Actualizar `odatasets_unique` con índices válidos
odatasets_unique <- odatasets_unique[!sapply(class_distributions, function(x) is.null(x$data)), ]
odatasets_unique$class_distribution <- sapply(valid_datasets, `[[`, "text")
odatasets_unique$class_distribution_vector <- lapply(valid_datasets, `[[`, "vector")
odatasets_unique$dataset <- lapply(valid_datasets, `[[`, "data")



#Datasets con solo variables independietes numericas

num_datasets <- odatasets_unique %>%
  filter(NumberOfSymbolicFeatures <= 1)# Filtrar datasets con NumberOfSymbolicFeatures <= 1


# Crear un data frame global para almacenar los resultados
global_results <- data.frame(
  name = character(),
  Best_Seed = integer(),
  ARI = numeric(),
  AMI = numeric(),
  NMI = numeric(),
  Mean_Silhouette = numeric(),
  Dunn_index = numeric(),
  Clusters = integer(),
  number_features = integer(),
  number_instances= integer(),
  cardinality_BAT = I(list()),
  cardinality_REAL = I(list()),
  stringsAsFactors = FALSE
)

# Registrar el tiempo inicial
start_time_total <- Sys.time()

# Iterar los datasets con manejo de errores
# Crear un vector vacío para almacenar los tiempos
execution_times <- numeric(nrow(odatasets_unique))
#nrow(odatasets_unique)
for (i in 1:nrow(odatasets_unique)) {
  cat("\n\n--- Ejecutando para dataset en posición:", i, "---\n")
  
  # Registrar el tiempo inicial
  start_time <- Sys.time()
  
  tryCatch({
    # Extraer el dataset
    dataset <- odatasets_unique[i]$dataset[[1]]
    dataset_name = odatasets_unique[i]$name
    # Verificar si el dataset es válido y preparar los datos
    prepared_data <- prepare_data(dataset)
    if (is.null(prepared_data)) {
      cat("Dataset en posición", i, "no tiene formato válido. Se omite.\n")
      execution_times[i] <- NA  # Guardar NA si se omite el dataset
      next
    }
    
    # Separar X e y
    X <- prepared_data$X
    y <- prepared_data$y
    
    # Extraer la cardinalidad objetivo
    target_cardinality <- odatasets_unique[i]$class_distribution_vector[[1]]
    
    # Verificar si la cardinalidad objetivo existe
    if (is.null(target_cardinality)) {
      cat("Cardinalidad objetivo no disponible para la posición", i, ". Se omite.\n")
      execution_times[i] <- NA  # Guardar NA si se omite el dataset
      next
    }
    
    # Ejecutar clustering
    run_clustering(dataset, target_cardinality, dataset_name)
    # Registrar el tiempo final
    end_time <- Sys.time()
    
    # Calcular la diferencia
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    execution_times[i] <- execution_time
    
    # Imprimir el tiempo de ejecución
    print(paste("Tiempo de ejecución para la posición", i, ":", execution_time, "segundos"))
    
  }, error = function(e) {
    # Manejar cualquier error y continuar con el siguiente dataset
    cat("Error al procesar el dataset en posición", i, ":", e$message, "\n")
    execution_times[i] <- NA  # Guardar NA en caso de error
  })
  
  # Registrar el tiempo final
  end_time <- Sys.time()
  
}


# Registrar el tiempo final
end_time_total <- Sys.time()

# Calcular la diferencia
execution_time_total <- end_time_total - start_time_total
print("Tiempo Total:")
print(execution_time_total)


execution_times <- execution_times[execution_times != 0 & !is.na(execution_times)]

#Resumen de Resultados

# Convertir listas en columnas para graficar
global_results$Cardinality_BAT <- sapply(global_results$cardinality_BAT, paste, collapse = ", ")
global_results$Cardinality_REAL <- sapply(global_results$cardinality_REAL, paste, collapse = ", ")

write.csv(global_results, "resultados.csv", row.names = FALSE)



# Preparar datos
violations_data <- data.frame(
  #Dataset = 1:nrow(global_results),
  Violations = sapply(1:nrow(global_results), function(i) {
    real <- unlist(global_results$cardinality_REAL[i])
    bat <- unlist(global_results$cardinality_BAT[i])
    sum(abs(real - bat))
  })
)

global_results_total <-  cbind(violations_data, global_results)

global_results_total$cardinality_BAT <- NULL
global_results_total$cardinality_REAL <- NULL

write.csv(global_results_total, "resultados.csv", row.names = FALSE)

