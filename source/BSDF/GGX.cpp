#include "GGX.h"

#define PI 3.1415926

// BRDF = F * D(H) * G(V, L) / (4.0f * dot(N, L) * dot(N, V))

bool GGX::SampleAnisoGGX(const glm::vec2& sample, const glm::vec3& viewDir,
                         glm::vec3& outDir, float& eval, float& pdf) {
  glm::vec3 localV = m_WorldToLocal * viewDir;

  glm::vec3 localH = SampleAnisoGGXNormal(sample, localV);
  glm::vec3 localL = glm::reflect(localV, localH);
  outDir = m_LocalToWorld * localL;

  if (localL.z < 0.001f || glm::dot(m_WorldNormal, outDir) <= 0.f) return false;

  float fresnel = 1.0f;

  eval = (fresnel * DAnisoGGX(localH) * SmithAnisoGGXShadowG2(localL, localV)) /
         (4.0f * localL.z * localV.z);

  pdf = (DAnisoGGX(localH) * abs(localH.z));

  eval *= localL.z;

  if (pdf < 0.001f) return false;
  return true;
}

// VNDF
// D_V(N) = G1(V)max(0, V * N)D(N) / (V * Z)
// PDF(L) = D_V(N) / (4 * (V * N))

bool GGX::SampleAnisoGGXVisible(const glm::vec2& sample,
                                const glm::vec3& viewDir, glm::vec3& outDir,
                                float& eval, float& pdf) {
  glm::vec3 localV = m_WorldToLocal * viewDir;

  glm::vec3 localH = SampleAnisoGGXVisibleNormal(sample, localV);

  glm::vec3 localL = glm::reflect(localV, localH);
  outDir = m_LocalToWorld * localL;

  if (localL.z < 0.001f || glm::dot(m_WorldNormal, outDir) <= 0.f) return false;

  float fresnel = 1.0f;

  eval = (fresnel * DAnisoGGX(localH) * SmithAnisoGGXShadowG2(localL, localV)) /
         (4.0f * localL.z * localV.z);

  pdf = (SmithAnisoGGXShadowG1(localV) * DAnisoGGX(localH)) / (4.0f * localV.z);

  eval *= localL.z;

  if (pdf < 0.001f) return false;
  return true;
}

bool GGX::SampleAnisoGGXVisibleSphericalCaps(const glm::vec2& sample,
                                             const glm::vec3& viewDir,
                                             glm::vec3& outDir, float& eval,
                                             float& pdf) {
  glm::vec3 localV = m_WorldToLocal * viewDir;

  glm::vec3 localH = SampleAnisoGGXVisibleNormalSphericalCaps(sample, localV);

  glm::vec3 localL = glm::reflect(localV, localH);
  outDir = m_LocalToWorld * localL;

  if (localL.z < 0.001f || glm::dot(m_WorldNormal, outDir) <= 0.f) return false;

  float fresnel = 1.0f;

  eval = (fresnel * DAnisoGGX(localH) * SmithAnisoGGXShadowG2(localL, localV)) /
         (4.0f * localL.z * localV.z);

  pdf = (SmithAnisoGGXShadowG1(localV) * DAnisoGGX(localH)) / (4.0f * localV.z);

  eval *= localL.z;

  if (pdf < 0.001f) return false;
  return true;
}

bool GGX::SampleAnisoGGXVisibleSphericalCapsBounded(const glm::vec2& sample,
                                                    const glm::vec3& viewDir,
                                                    glm::vec3& outDir,
                                                    float& eval, float& pdf) {
  glm::vec3 localV = m_WorldToLocal * viewDir;

  glm::vec3 localH =
      SampleAnisoGGXVisibleNormalSphericalCapsBounded(sample, localV);

  glm::vec3 localL = glm::reflect(localV, localH);
  outDir = m_LocalToWorld * localL;

  if (localL.z < 0.001f || glm::dot(m_WorldNormal, outDir) <= 0.f) return false;

  float fresnel = 1.0f;

  eval = (fresnel * DAnisoGGX(localH) * SmithAnisoGGXShadowG2(localL, localV)) /
         (4.0f * localL.z * localV.z);

  pdf = GGXBoundedVndfPDF(localL, localV);

  eval *= localL.z;

  if (pdf < 0.001f) return false;
  return true;
}

glm::vec3 GGX::SampleAnisoGGXNormal(const glm::vec2& sample,
                                    const glm::vec3& localV) {
  float cosTheta = 0.0f;
  float phi = 2 * PI * sample.y;

  if (m_Roughness.x == m_Roughness.y) {
    float tanTheta2 =
        m_Roughness.x * m_Roughness.x * sample.x / (1.0f - sample.x);
    cosTheta = 1.0f / sqrtf(1.0f + tanTheta2);
  } else {
    phi = atanf(m_Roughness.y / m_Roughness.x *
                tanf(2 * PI * sample.y + 0.5f * PI));
    if (sample.y > 0.5f) phi += PI;
    float sinPhi = sinf(phi);
    float cosPhi = cosf(phi);
    const float alphax2 = m_Roughness.x * m_Roughness.x;
    const float alphay2 = m_Roughness.y * m_Roughness.y;
    const float alpha2 =
        1.0f / (cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
    float tanTheta2 = alpha2 * sample.x / (1.0f - sample.x);
    cosTheta = 1.0f / sqrtf(1.0f + tanTheta2);
  }
  float sinTheta = sqrtf(glm::max(0.0f, 1.0f - cosTheta * cosTheta));
  glm::vec3 Ne = glm::normalize(
      glm::vec3(sinTheta * cosf(phi), sinTheta * sinf(phi), cosTheta));
  return Ne;
}

glm::vec3 GGX::SampleAnisoGGXVisibleNormal(const glm::vec2& sample,
                                           const glm::vec3& localV) {
  // transforming the view direction to the hemisphere configuration
  glm::vec3 Vh = glm::normalize(
      glm::vec3(m_Roughness.x * localV.x, m_Roughness.y * localV.y, localV.z));
  // orthonormal basis(with special case if cross product is zero)
  float lensq = Vh.x * Vh.x + Vh.y * Vh.y;
  glm::vec3 T1 = lensq > 0.0f ? glm::vec3(-Vh.y, Vh.x, 0.0f) / sqrtf(lensq)
                              : glm::vec3(1.0f, 0.0f, 0.0f);
  glm::vec3 T2 = glm::cross(Vh, T1);

  // parameterization of the projected area
  float r = sqrtf(sample.x);
  float phi = 2.0f * PI * sample.y;
  float t1 = r * cos(phi);
  float t2 = r * sin(phi);
  float s = 0.5 * (1.0f + Vh.z);
  t2 = (1.0f - s) * sqrtf(1.0 - t1 * t1) + s * t2;

  // reprojection onto hemisphere
  glm::vec3 Nh =
      t1 * T1 + t2 * T2 + sqrtf(glm::max(0.0f, 1.0f - t1 * t1 - t2 * t2)) * Vh;

  // transforming the normal back to the ellipsoid configuration
  glm::vec3 Ne = glm::normalize(glm::vec3(
      m_Roughness.x * Nh.x, m_Roughness.y * Nh.y, glm::max(0.0f, Nh.z)));

  return Ne;
}

glm::vec3 GGX::SampleAnisoGGXVisibleNormalSphericalCaps(
    const glm::vec2& sample, const glm::vec3& localV) {
  // transforming the view direction to the hemisphere configuration
  glm::vec3 Vh = glm::normalize(
      glm::vec3(m_Roughness.x * localV.x, m_Roughness.y * localV.y, localV.z));

  // sample a spherical cap in (-Vh.z 1]
  float phi = 2.0f * PI * sample.x;
  // float z = (1.0f - sample.y) * (1.0f + localV.z) - localV.z;
  float z = 1.0 - sample.y - localV.z * sample.y;
  float sinTheta = sqrtf(glm::clamp(1.0f - z * z, 0.0f, 1.0f));
  float x = sinTheta * cos(phi);
  float y = sinTheta * sin(phi);
  glm::vec3 c = glm::vec3(x, y, z);

  // compute halfway direction
  glm::vec3 h = c + Vh;

  // transforming the normal back to the ellipsoid configuration
  glm::vec3 Ne = glm::normalize(
      glm::vec3(m_Roughness.x * h.x, m_Roughness.y * h.y, glm::max(0.0f, h.z)));

  return Ne;
}

glm::vec3 GGX::SampleAnisoGGXVisibleNormalSphericalCapsBounded(
    const glm::vec2& sample, const glm::vec3& localV) {
  // transforming the view direction to the hemisphere configuration
  glm::vec3 Vh = glm::normalize(
      glm::vec3(m_Roughness.x * localV.x, m_Roughness.y * localV.y, localV.z));

  // Sample a spherical cap
  float phi = 2.0f * PI * sample.x;
  float a = glm::clamp(glm::min(m_Roughness.x, m_Roughness.y), 0.0f, 1.0f);
  float s = 1.0f + glm::length(glm::vec2(localV.x, localV.y));
  float a2 = a * a;
  float s2 = s * s;
  float k = (1.0f - a2) * s2 / (s2 + a2 * localV.z * localV.z);
  float b = localV.z > 0.0f ? k * Vh.z : Vh.z;
  // float z = (1.0f - sample.y) * (1.0f + b) - b;
  float z = 1.0f - sample.y - b * sample.y;
  float sinTheta = sqrtf(glm::clamp(1.0f - z * z, 0.0f, 1.0f));
  float x = sinTheta * cos(phi);
  float y = sinTheta * sin(phi);
  glm::vec3 c = glm::vec3(x, y, z);

  // compute halfway direction
  glm::vec3 h = c + Vh;

  // transforming the normal back to the ellipsoid configuration
  glm::vec3 Ne = glm::normalize(
      glm::vec3(m_Roughness.x * h.x, m_Roughness.y * h.y, glm::max(0.0f, h.z)));

  return Ne;
}

float GGX::DAnisoGGX(const glm::vec3& normal) {
  float x = normal.x / m_Roughness.x;
  float y = normal.y / m_Roughness.y;
  float z = normal.z;

  float alpha = x * x + y * y + z * z;

  return 1.0f / (PI * m_Roughness.x * m_Roughness.y * alpha * alpha);
}

float GGX::SmithAnisoGGXShadowG1(const glm::vec3& viewDir) {
  return 1.0f / (1.0f + Lambda(viewDir));
}

float GGX::SmithAnisoGGXShadowG2(const glm::vec3& wi, const glm::vec3& wo) {
  return 1.0f / (1.0f + Lambda(wi) + Lambda(wo));
}

float GGX::Lambda(const glm::vec3& viewDir) {
  float x = m_Roughness.x * viewDir.x;
  float y = m_Roughness.y * viewDir.y;
  float z = viewDir.z;

  float alpha = (x * x + y * y) / (z * z);
  return (-1.0f + sqrtf(1.0f + alpha)) / 2.0f;
}

float GGX::GGXBoundedVndfPDF(const glm::vec3& wi, const glm::vec3& wo) {
  glm::vec3 wm = glm::normalize(wi + wo);
  float ndf = DAnisoGGX(wm);
  glm::vec2 ai = glm::vec2(m_Roughness.x * wi.x, m_Roughness.y * wi.y);
  float len2 = glm::dot(ai, ai);
  glm::vec3 Vh = glm::normalize(glm::vec3(ai, wi.z));
  float t = glm::length(Vh);
  if (wi.z >= 0.0f) {
    float alpha =
        glm::clamp(glm::min(m_Roughness.x, m_Roughness.y), 0.0f, 1.0f);
    float s = 1.0f + glm::length(glm::vec2(wi.x, wi.y));
    float alpha2 = alpha * alpha;
    float s2 = s * s;
    float k = (1.0f - alpha2) * s2 / (s2 + alpha2 * wi.z * wi.z);
    return ndf / (2.0f * (k * wi.z + t));
  }
  return ndf * (t - wi.z) / (2.0f * len2);
}