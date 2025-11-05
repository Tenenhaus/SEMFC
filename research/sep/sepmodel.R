
#lavaan

sep.model <-  '
# latent variable definitions
clinique =~ time+errors+front+all+speed
biologique =~ polig+micro+mog+ftl+pmog+pftl


# Regressions
clinique ~ biologique

'


#semfc

C <- matrix(c(0, 0,
              1, 0), 2, 2, byrow = TRUE)

colnames(C) <- rownames(C) <- c('clinique', 'biologique')

mode <- rep("reflective", 2)

