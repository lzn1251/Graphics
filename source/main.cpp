

#include <iostream>

#include "BSDF/GGX.h"

#define N 99

int main() {
  glm::vec2 sample;
  glm::vec3 localV = glm::normalize(glm::vec3(0.8f, 0.0f, 0.2f));
  GGX ggx(glm::vec2(0.8f, 0.8f));
  srand(time(NULL));
  unsigned int ggxHCount = 0;
  unsigned int ggxVisHCount = 0;
  unsigned int ggxScHCount = 0;
  unsigned int ggxBoundedCount = 0;
  for (int i = 0; i < 10000; ++i) {
    sample.x = rand() % (N + 1) / (float)(N + 1);
    sample.y = rand() % (N + 1) / (float)(N + 1);
    // std::cout << sample.x << " " << sample.y << std::endl;

    glm::vec3 ggxH = ggx.SampleAnisoGGXNormal(sample, localV);
    if (glm::dot(ggxH, localV) > 0.0f) ggxHCount++;

    glm::vec3 ggxVisH = ggx.SampleAnisoGGXVisibleNormal(sample, localV);
    if (glm::dot(ggxVisH, localV) > 0.0f) ggxVisHCount++;

    glm::vec3 ggxScH =
        ggx.SampleAnisoGGXVisibleNormalSphericalCaps(sample, localV);
    if (glm::dot(ggxScH, localV) > 0.0f) ggxScHCount++;

    glm::vec3 ggxBoundedH =
        ggx.SampleAnisoGGXVisibleNormalSphericalCapsBounded(sample, localV);
    if (glm::dot(ggxBoundedH, localV) > 0.0f) ggxBoundedCount++;
  }

  std::cout << std::endl;

  std::cout << ggxHCount << std::endl;
  std::cout << ggxVisHCount << std::endl;
  std::cout << ggxScHCount << std::endl;
  std::cout << ggxBoundedCount << std::endl;
}